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

    # interface flux
    uf = rd.Vf * u
    uP = uf[md.mapP]

    interface_flux = @. surface_flux(uf, uP, SVector{1}(md.nxJ), equations)
    du .= rd.LIFT * interface_flux
   
    rhs_vol_high = similar(du, size(du, 1))
    rhs_vol_low = similar(rhs_vol_high)
    rhs_vol_hyper = similar(rhs_vol_high)
    
    A = fill!(similar(du, size(du, 1), size(du, 1)), zero(eltype(du)))
    a = zeros(size(du, 1), size(du, 1))
    for e in 1:md.num_elements        
        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))
        fill!(rhs_vol_hyper, zero(eltype(rhs_vol_hyper)))
        
        # high order update
        u_element = view(u, :, e)
        v = cons2entropy.(u_element, equations)
        for i in eachindex(u_element), j in eachindex(u_element)            
            if i > j
                fij = volume_flux(u_element[i], u_element[j], 1, equations)

                FH_ij = Q_skew[i, j] * fij
                rhs_vol_high[i] +=  FH_ij
                rhs_vol_high[j] += -FH_ij
                
                # True low order stuffs
                if blend != :ELXF
                    FL_ij = zero(eltype(rhs_vol_low))
                    if abs(Q_skew_low[i, j]) > 1e-12
                        nij = Q_skew_low[i, j] / abs(Q_skew_low[i, j])
                        fij = surface_flux(u_element[i], u_element[j], SVector{1}(nij), equations)
                        FL_ij = abs(Q_skew_low[i, j]) * fij
                        rhs_vol_low[i] += FL_ij
                        rhs_vol_low[j] -= FL_ij    
                    end
                end

                # Fake low order stuffs
                if blend == :ELXF || blend == :NLXF
                    FL_ij = zero(eltype(rhs_vol_low))
                    if abs(Q_skew_low[i, j]) > 1e-12
                        fij = u_element[i] - u_element[j]
                        FL_ij = abs(Q_skew_low[i, j]) * fij
                        rhs_vol_low[i] += FL_ij
                        rhs_vol_low[j] -= FL_ij    
                    end
                end
                if blend == :NLXF
                    FL_ij = zero(eltype(rhs_vol_low))
                    if abs(Q_skew_hyper[i, j]) > 1e-12
                        fij = u_element[i] - u_element[j]
                        FL_ij = abs(Q_skew_hyper[i, j]) * fij
                        rhs_vol_hyper[i] += FL_ij
                        rhs_vol_hyper[j] -= FL_ij
                    end
                end
            end
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

        if blend == :ELXF
            l = sum(dot.(v, rhs_vol_low))
            dΨ = sum(cache.B * psi.(u_element, equations))
            Ec = sum(dot.(v, rhs_vol_high))
            λ = 0

            if abs(l) > 100 * eps() && (-Ec > dΨ + 100 * eps())
                λ = 2(dΨ + Ec) / l
            end

            r = (rhs_vol_high - λ / 2 * rhs_vol_low)

            # if !(-sum(dot.(v, r)) <= dΨ + 300 * eps())
            #     println("VIOLATION: ", -sum(dot.(v, r)) - dΨ)
            # end

            view(du, :, e) .+= rd.M \ r
        elseif blend == :NLXF
            l = sum(dot.(v, rhs_vol_low))
            l_hyper = sum(dot.(v, rhs_vol_hyper))
            dΨ = sum(cache.B * psi.(u_element, equations))
            Ec = sum(dot.(v, rhs_vol_high))
            λ = 0
            b = dΨ + Ec

            a = [l, l_hyper]

            # @assert l > -1000 * eps()
            # @assert l_hyper > -100000 * eps()

            if !(l_hyper > -100 * eps()) && !isnan(l_hyper)
                println(l_hyper)
            end

            # Need aTθ = b

            θ = [0, 0]

            if norm(a) > 100 * eps() && (-Ec > dΨ + 100 * eps())
                θ = b * a / a'a
            end

            # println(θ)

            r = (rhs_vol_high - θ[1] * rhs_vol_low - θ[2] * rhs_vol_hyper)

            view(du, :, e) .+= rd.M \ r
        elseif blend == :elementwise
            # 1. elementwise blending of rhs_vol_high/low
            denom = -(sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(1.0, θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)
        elseif blend == :subcell
            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))

            # Call the Knapsack Solver
            a = vec(a)
            θ = cache.knapsack_solver(a, b)
            
            if nodewise_shock_capturing > 0
                @. a *= nodewise_shock_capturing * (θ - 1)^2 + 1

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
                @. a *= nodewise_shock_capturing * (θ - 1)^2 + 1

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
    
    # invert Jacobian and mass matrix
    du ./= -md.J
end

rd = RefElemData(Line(), SBP(), N)
(VX, ), EToV = uniform_mesh(Line(), K)

md = MeshData((VX,), EToV, rd)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

# subcell FV operators
Δ = diagm(0 => -ones(rd.N+1), 1=>ones(rd.N+1))[1:end-1,:]
R = tril(ones(size(Δ')), -1)