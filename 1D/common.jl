import Trixi: flux_hllc
# use the 2D implementation since the 1D version doesn't account for n::Integer < 0
function flux_hllc(u_ll, u_rr, n::SVector{1}, equations::CompressibleEulerEquations1D)
    f = flux_hllc(SVector(u_ll[1], u_ll[2], 0, u_ll[3]), 
                  SVector(u_rr[1], u_rr[2], 0, u_rr[3]), 
                  SVector(n[1], 0.0), 
                  CompressibleEulerEquations2D(equations.gamma))    
    return SVector(f[1], f[2], f[4])
end

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

        function C(u)
            if u >= 0
                return min(u, one(typeof(u)))
            else
                return one(typeof(u))
            end
        end

        function preserve_positivity_backward_nodewise(type, component)
            if preserve_positivity >= 0
                # First, determine h
                h = (1 - preserve_positivity) * (1 / current_timestep * md.J[1, 1] * rd.M * u_element - rd.M * view(du, :, e))
                # Now, \hat{h}
                h_hat = getindex.(h - (1 - preserve_positivity) * rhs_vol_low, component)
                # Now, c
                c = getindex.(R * (rhs_vol_high - rhs_vol_low), component)

                # Now, we determine 1 - l_c
                l_c = Vector{type}(undef, size(R, 1))
                l_c[1] = C(h_hat[1] / (-2c[1]))
                l_c[N + 2] = C(h_hat[end] / (2c[end]))
                for i in 2:(size(R, 1) - 1)
                    l_c[i] = min(C(h_hat[i - 1] / (2c[i])), C(h_hat[i] / (-2c[i])))
                end

                return 1 .- l_c
            else
                return zeros(size(R, 1))
            end
        end

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
            maximum_theta = preserve_positivity_forward()

            denom = -(sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(maximum_theta, θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)
        elseif blend == :subcell
            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))

            # Call the Knapsack Solver
            a = vec(a)
            θ = cache.knapsack_solver(a, b, upper_bounds = ones(eltype(a), length(a))) # upper_bounds = exp.(.005 * [0, 0, -1, -2, -3, -4, -5, -6]) # upper_bounds = ones(eltype(a), length(a)) exp.(.002 * [0, 0, -1, -2])
            
            if nodewise_shock_capturing > 0
                @. a *= N * tan(pi/2 * nodewise_shock_capturing)^2 * (1 - θ)^2 + 1

                θ = cache.knapsack_solver(a, b)
            end

            view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ) * R * rhs_vol_high + Diagonal(1 .- θ) * R * rhs_vol_low))
        elseif blend == :subcell_reversed
            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))

            # minimum_theta = preserve_positivity_backward_again()
            l_c_1 = preserve_positivity_backward_nodewise(eltype(a), 1)
            l_c_3 = preserve_positivity_backward_nodewise(eltype(a), 3)
            l_c = max.(l_c_1, l_c_3)
            if initial_condition == initial_condition_leblanc_shocktube && t < 1e-7
                println("Elementwise limiting coefficients chosen for lablanc")
                l_c .= maximum(l_c)
                l_c = min.(1, 2l_c)
            end

            # l_c = 1 .- exp.(-.02 * [0, (0:(size(R, 1) - 2))...])
            # Call the Knapsack Solver
            a = vec(a)
            b = -sum(cache.B * psi.(u_element, equations)) - sum(dot.(v, rhs_vol_high)) - a'l_c # minimum_theta * sum(a)
            θ = cache.knapsack_solver(a, b, upper_bounds = (1 .- l_c)) # ones(length(a)) * (1 - minimum_theta)
            
            if nodewise_shock_capturing > 0
                @. θ *= N^(nodewise_shock_capturing) * θ^2 / ((1 - θ)^2 + 1e-10) + 1
                θ = clamp.(θ, 0, 1 .- l_c)
            end

            view(du, :, e) .+= rd.M \ (rhs_vol_high + Δ * Diagonal(θ + l_c) * R * (rhs_vol_low - rhs_vol_high)) # θ .+ minimum_theta
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
    
    # invert Jacobian and mass matrix
    if any(isnan.(sum(du)))
        println("NAN DETECTED")
    end
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

function niceplot!()
    plot!(legendfontsize=12)
    plot!(xtickfontsize=12)
    plot!(ytickfontsize=12)
    plot!(labelfontsize=14)
end

current_timestep = dt # for positivity preservation
if knapsack_stats != Nothing # if we are grabbing stats, reset them on new run
    knapsack_stats = []
end