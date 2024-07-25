function rhs!(du, u, cache, t)
    (; physics, rd, md) = cache
    (; Q_skew) = cache.high_order_operators
    (; Q_skew_low) = cache.low_order_operators    
    (; Δ, R) = cache.fv_operators
    (; equations, volume_flux, surface_flux) = physics
    (; bc, VDM_inv, shock_capturing, nodewise_shock_capturing) = cache

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
    
    A = fill!(similar(du, size(du, 1), size(du, 1)), zero(eltype(du)))
    a = zeros(size(du, 1), size(du, 1))
    for e in 1:md.num_elements        
        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))
        
        # high order update
        u_element = view(u, :, e)
        v = cons2entropy.(u_element, equations)
        for i in eachindex(u_element), j in eachindex(u_element)            
            if i > j
                fij = volume_flux(u_element[i], u_element[j], 1, equations)

                FH_ij = Q_skew[i, j] * fij
                rhs_vol_high[i] +=  FH_ij
                rhs_vol_high[j] += -FH_ij

                FL_ij = zero(eltype(rhs_vol_low))
                if abs(Q_skew_low[i, j]) > 1e-12
                    nij = Q_skew_low[i, j] / abs(Q_skew_low[i, j])
                    fij = surface_flux(u_element[i], u_element[j], SVector{1}(nij), equations)
                    FL_ij = abs(Q_skew_low[i, j]) * fij
                    rhs_vol_low[i] += FL_ij
                    rhs_vol_low[j] -= FL_ij    
                end

                # create antidiffusive flux matrix
                A[i,j] = FH_ij - FL_ij
                A[j,i] = -A[i, j]
            end
        end
        
        # optimization target
        b = sum(cache.B * psi.(u_element, equations)) + sum(dot.(v, rhs_vol_low))

        ### Shock capturing optimization target
        if shock_capturing
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

            b = (1 - .01 * ϵ) * b
        end

        (; blend) = cache

        if blend == :elementwise
            # 1. elementwise blending of rhs_vol_high/low
            denom = (sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(1.0, θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)
        elseif blend == :subcell
            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))

            # Call the Knapsack Solver
            a = vec(a)
            θ = cache.knapsack_solver(a, b)
            
            if nodewise_shock_capturing
                epsilon = 240.0
                # @. a *= -log(1/e * θ)
                # @. a *= 1 / cbrt(θ)
                @. a *= 1 - 4 * N^2 * log(min(θ + 1 / (2N * K), 1))
                # # @. a *= (1 - 2 * N * (K / 64) * log(θ + 1e-8))
                # @. a *= -epsilon * θ + (1 + epsilon)
                # # @. a *= (-sqrt(epsilon) * θ + sqrt(epsilon))^2 + 1
                # @. a *= -epsilon * min(1, θ + 1e-3) + (1 + epsilon)

                θ = cache.knapsack_solver(a, b)
                # @. a *= -epsilon * θ + (1 + epsilon)
                # # @. a *= (-sqrt(epsilon) * θ + sqrt(epsilon))^2 + 1

                # θ = cache.knapsack_solver(a, b)
            end
            view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ) * R * rhs_vol_high + Diagonal(1 .- θ) * R * rhs_vol_low))
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