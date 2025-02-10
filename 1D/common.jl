function graph_laplacian!(vec, matr, operator; raw=false)
    for i in eachindex(vec), j in eachindex(vec)
        if i > j
            FL_ij = zero(eltype(vec))
            if abs(matr[i, j]) > 1e-12
                fij = operator(i, j)
                comp = matr[i, j]
                if !raw
                    comp = abs(comp)
                end
                FL_ij = comp * fij
                vec[i] += FL_ij
                vec[j] -= FL_ij
            end
        end
    end
end


function rhs!(du, u, cache, t)
    (; physics, rd, md) = cache
    (; Q_skew) = cache.high_order_operators
    (; Q_skew_low) = cache.low_order_operators
    (; Q_skew_hyper) = cache.hyper_operator
    (; Δ, R) = cache.fv_operators
    (; equations, volume_flux, surface_flux) = physics
    (; bc, VDM_inv, shock_capturing, nodewise_shock_capturing) = cache
    (; blend) = cache

    # Boundary conditions (set end points to initial conditions)
    if !isnothing(bc)
        u[1, 1] = bc[:, 1]
        u[end, end] = bc[:, 2]
    end

    viscosities = Vector{Float64}(undef, K)

    # interface flux
    uf = rd.Vf * u
    uP = uf[md.mapP]

    interface_flux = @. surface_flux(uf, uP, SVector{1}(md.nxJ), equations)
    du .= rd.LIFT * interface_flux
   
    rhs_vol_high = similar(du, size(du, 1))
    rhs_vol_low = similar(rhs_vol_high)
    rhs_vol_visc = similar(rhs_vol_high)
    
    A = fill!(similar(du, size(du, 1), size(du, 1)), zero(eltype(du)))
    for e in 1:md.num_elements        
        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))
        fill!(rhs_vol_visc, zero(eltype(rhs_vol_visc)))
        
        # high order update
        u_element = view(u, :, e)
        v = cons2entropy.(u_element, equations)

        high_order_operator(i, j) = volume_flux(u_element[i], u_element[j], 1, equations)
        function low_order_operator(i, j)
            nij = Q_skew_low[i, j] / abs(Q_skew_low[i, j])
            return surface_flux(u_element[i], u_element[j], SVector{1}(nij), equations)
        end

        graph_laplacian!(rhs_vol_high, Q_skew, high_order_operator; raw=true)
        
        if blend != :viscosity && blend != :nodewiseviscosity # the only methods which do not need rhs_vol_low are viscosity and nodewise viscosity
            graph_laplacian!(rhs_vol_low, Q_skew_low, low_order_operator, raw=false)
        end

        # The only methods which need a viscous term are the viscosity methods
        if blend == :viscosity || blend == :hyperviscosity || blend == :nodewiseviscosity
            low_order_operator_diffusive(i, j) = u_element[i] - u_element[j]
            graph_laplacian!(rhs_vol_visc, Q_skew_low, low_order_operator_diffusive, raw=false)
        end
        
        # optimization target
        b = sum(cache.B * psi.(u_element, equations)) + sum(dot.(v, rhs_vol_low))

        ### Shock capturing optimization target
        if shock_capturing > 0
            μ = VDM_inv * getindex.(u_element, 1)
            sum_to_end = sum(μ.^2)

            s = log10(max(μ[end]^2 / sum_to_end, μ[end - 1]^2 / (sum_to_end - μ[end]^2)))
            κ = 1.0
            s0 = -4 * log10(size(u_element, 1) - 1)
            ϵ = 0.0
            if s < s0 + κ
                ϵ = 1.0
            elseif s < s0 - κ
                ϵ = 0.0
            else
                ϵ = 0.5 * (1.0 - sin(π * (s - s0) / (2 * κ)))
            end

            b = (1 - shock_capturing * ϵ) * b
        end

        if blend == :viscosity
            l = sum(dot.(v, rhs_vol_visc))
            dΨ = sum(cache.B * psi.(u_element, equations))
            Ec = sum(dot.(v, rhs_vol_high))
            λ = 0

            if abs(l) > 100 * eps() && (-Ec > dΨ + 100 * eps())
                λ = 2(dΨ + Ec) / l
            end

            r = (rhs_vol_high - λ / 2 * rhs_vol_visc)

            viscosities[e] = λ / 2

            # if !(-sum(dot.(v, r)) <= dΨ + 300 * eps())
            #     println("VIOLATION: ", -sum(dot.(v, r)) - dΨ)
            # end

            view(du, :, e) .+= rd.M \ r
        elseif blend == :hyperviscosity
            w0 = 1
            w1 = 1
            l = sum(dot.(v, rhs_vol_visc))
            l_yimin = sum(dot.(v, rhs_vol_low))
            dΨ = sum(cache.B * psi.(u_element, equations))
            Ec = sum(dot.(v, rhs_vol_high))
            λ = 0
            b = dΨ + Ec

            a = [w0 * l, w1 * l_yimin]

            θ = [0, 0, 0]

            if norm(a) > 100 * eps() && (-Ec > dΨ + 100 * eps())
                θ = b * a / a'a
            end

            r = (rhs_vol_high - w0 * θ[1] * rhs_vol_visc - w1 * θ[2] * rhs_vol_low)

            view(du, :, e) .+= rd.M \ r
        elseif blend == :nodewiseviscosity
            a = vec(v' * Δ * Diagonal(R * rhs_vol_visc))
            ac = clamp.(a, 0, Inf)

            dΨ = sum(cache.B * psi.(u_element, equations))
            Ec = sum(dot.(v, rhs_vol_high))

            λ = 0 * a
            if a'ac > 100 * eps() && (-Ec > dΨ + 100 * eps())
                b = dΨ + Ec
                λ = b * ac / (a'ac)
            end

            r = rhs_vol_high - Δ * (Diagonal(λ) * R * rhs_vol_visc)

            view(du, :, e) .+= rd.M \ r
        elseif blend == :elementwise
            # 1. elementwise blending of rhs_vol_high/low
            denom = -(sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(1.0, θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)
        elseif blend == :subcell
            # Here we change the entropy variables and psi
            # v .= u_element
            
            # println(flux.(u_element, 1, equations)' * rd.M * rd.Dr * u_element)
            # dPsi = sum(flux.(u_element, 1, equations)' * rd.M * rd.Dr * u_element)
            # b = dPsi + sum(dot.(v, rhs_vol_low))


            # Here is the end of that



            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))

            # Call the Knapsack Solver
            a = vec(a)
            θ = cache.knapsack_solver(a, b)
            
            if nodewise_shock_capturing > 0
                @. a *= N * tan(pi/2 * nodewise_shock_capturing)^2 * (1 - θ)^2 + 1

                θ = cache.knapsack_solver(a, b)
            end
            view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ) * R * rhs_vol_high + Diagonal(1 .- θ) * R * rhs_vol_low))
        elseif blend == :subcell_dynamic
            a = map((a, b) -> a .* b, Δ' * v, (R * (rhs_vol_low - rhs_vol_high)))
            a = reinterpret(Float64, a)

            # Call the Knapsack Solver
            a = vec(a)
            θ = cache.knapsack_solver(a, b)

            if nodewise_shock_capturing > 0
                @. a *= N * tan(pi/2 * nodewise_shock_capturing)^2 * (1 - θ)^2 + 1

                θ = cache.knapsack_solver(a, b)
            end

            θ = reinterpret(SVector{3, Float64}, θ)
            
            view(du, :, e) .+= rd.M \ (Δ * (map((θ, a, b) -> θ .* a + (SVector{3, Float64}(1, 1, 1) - θ) .* b, θ, R * rhs_vol_high, R * rhs_vol_low)))
        elseif blend == :loworder
            view(du, :, e) .+= rd.M \ rhs_vol_low
        else
            view(du, :, e) .+= rd.M \ rhs_vol_high
        end
    end

    if blend == :viscosity
        # Add the surface term back
        interface_flux = ((uf - uP) .* md.nxJ) * Diagonal(viscosities)
        du .-= rd.LIFT * interface_flux
    end
    
    # invert Jacobian and mass matrix
    du ./= -md.J

    # if any(isnan.(sum.(du)))
    #     println("NAN DETECTED")
    # end
end

rd = RefElemData(Line(), SBP(), N)
(VX, ), EToV = uniform_mesh(Line(), K)

md = MeshData((VX,), EToV, rd)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

# subcell FV operators
Δ = diagm(0 => -ones(rd.N+1), 1=>ones(rd.N+1))[1:end-1,:]
R = tril(ones(size(Δ')), -1)