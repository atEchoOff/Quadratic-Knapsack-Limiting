function rhs!(du, u, cache, t)
    (; rd, md) = cache
    (; Qr_skew, Qs_skew) = cache
    (; Qr_sparse_skew, Qs_sparse_skew, Δr, Rr) = cache
    (; equations, volume_flux, surface_flux) = cache
    (; VDM_inv, shock_capturing) = cache
    (; rhs_vol_high_cache, rhs_vol_low_cache, rhs_vol_visc_cache, v_cache) = cache

    # interface flux
    uf = rd.Vf * u
    uP = uf[md.mapP]
    interface_flux = @. surface_flux(uf, uP, 
                                     SVector(md.nx, md.ny), equations) * md.Jf
    du .= rd.LIFT * interface_flux
   
    @inbounds Threads.@threads for e in axes(du, 2)
    # for e in axes(du, 2)
        rhs_vol_high = view(rhs_vol_high_cache, :, e)
        rhs_vol_low = view(rhs_vol_low_cache, :, e)
        rhs_vol_visc = view(rhs_vol_visc_cache, :, e)
        v = view(v_cache, :, e)

        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))
        fill!(rhs_vol_visc, zero(eltype(rhs_vol_visc)))

        # if blend == :viscosity
        #     # Add surface contribution to r, and zero out du
        #     rhs_vol_high .= view(du, :, e)
        #     fill!(view(du, :, e), zero(eltype(du)))
        # end

        u_element = view(u, :, e)
        v .= cons2entropy.(u_element, equations)
        # v_avg = cons2entropy(.5 * sum(rd.wq .* u_element), equations)

        for i in axes(du, 1), j in axes(du, 1)
            if i > j
                # high order
                Qx_ij = md.rxJ[1,e] * Qr_skew[i,j] + md.sxJ[1,e] * Qs_skew[i,j]
                Qy_ij = md.ryJ[1,e] * Qr_skew[i,j] + md.syJ[1,e] * Qs_skew[i,j]
                rhs_ij = zero(eltype(du))
                if abs(Qx_ij) > 100 * eps() 
                    Fx_ij = volume_flux(u[i,e], u[j,e], 1, equations)
                    rhs_ij += Qx_ij * Fx_ij
                end
                if abs(Qy_ij) > 100 * eps()                    
                    Fy_ij = volume_flux(u[i,e], u[j,e], 2, equations)
                    rhs_ij += Qy_ij * Fy_ij
                end
                rhs_vol_high[i] += rhs_ij
                rhs_vol_high[j] -= rhs_ij

                # No need for low order du for purely viscosity methods
                if blend != :viscosity && blend != :nodewiseviscosity
                    Qx_ij = md.rxJ[1,e] * Qr_sparse_skew[i,j] + md.sxJ[1,e] * Qs_sparse_skew[i,j]
                    Qy_ij = md.ryJ[1,e] * Qr_sparse_skew[i,j] + md.syJ[1,e] * Qs_sparse_skew[i,j]
                    nij_scaled = SVector(Qx_ij, Qy_ij)
                    nij_norm = norm(nij_scaled)
                    if nij_norm > 100 * eps()
                        nij = nij_scaled / nij_norm
                        Fij = surface_flux(u[i, e], u[j, e], nij, equations)
                        rhs_low_ij = Fij * nij_norm
                        rhs_vol_low[i] += rhs_low_ij
                        rhs_vol_low[j] -= rhs_low_ij
                    end
                end

                # Only need viscosity term for viscosity methods
                if blend == :viscosity || blend == :hyperviscosity || blend == :nodewiseviscosity
                    Qx_ij = md.rxJ[1,e] * Qr_sparse_skew[i,j] + md.sxJ[1,e] * Qs_sparse_skew[i,j]
                    Qy_ij = md.ryJ[1,e] * Qr_sparse_skew[i,j] + md.syJ[1,e] * Qs_sparse_skew[i,j]
                    nij_scaled = SVector(Qx_ij, Qy_ij)
                    nij_norm = norm(nij_scaled)
                    if nij_norm > 100 * eps()
                        ###FIXME derivative calculation
                        # dudv = (ForwardDiff.jacobian(x -> entropy2cons(x, equations), v_avg))

                        ###

                        Fij = u[i, e] - u[j, e]
                        # Fij = dudv * (v[i] - v[j])

                        rhs_low_ij = Fij * nij_norm
                        rhs_vol_visc[i] += rhs_low_ij
                        rhs_vol_visc[j] -= rhs_low_ij
                    end
                end
            end
        end

        function preserve_positivity_forward()
            maximum_theta = 1

            if preserve_positivity >= 0
                # We need to determine the bound here. 
                # top part first
                top_part = (1 - preserve_positivity) * (rd.M * (md.J[1, 1] / current_timestep * u_element - view(du, :, e)) - rhs_vol_low)

                # Now bottom part
                bottom_part = rhs_vol_high - rhs_vol_low

                # FIXME for compressible euler, we satisfy positivity just for 1st and 4th component...
                top_part = vcat(getindex.(top_part, 1), getindex.(top_part, 4))
                bottom_part = vcat(getindex.(bottom_part, 1), getindex.(bottom_part, 4))

                # If bottom part is negative, we dont want to count it in the minimum.
                top_part = top_part[bottom_part .> 0]
                bottom_part = 2 * bottom_part[bottom_part .> 0]

                # Now, compute the minimum
                if length(bottom_part) > 0
                    maximum_theta = minimum(top_part ./ bottom_part)
                    maximum_theta = min(1, maximum_theta)
                end
            end

            return maximum_theta
        end

        function C(u)
            if u >= 0
                return min(u, 1)
            else
                return 1
            end
        end

        function preserve_positivity_backward_nodewise(component)
            if preserve_positivity >= 0
                # First, determine h
                h = (1 - preserve_positivity) * (1 / current_timestep * md.J[1, 1] * rd.M * u_element - rd.M * view(du, :, e))
                # Now, \hat{h}
                h_hat = getindex.(h - (1 - preserve_positivity) * rhs_vol_low, component)
                # Now, c
                c = getindex.(Rr * (rhs_vol_high - rhs_vol_low), component)

                # Now, we determine 1 - l_c
                l_c = Vector{Float64}(undef, size(Rr, 1))
                l_c[1] = C(h_hat[1] / (-2c[1]))
                l_c[end] = C(h_hat[end] / (2c[end]))
                for i in 2:(size(Rr, 1) - 1)
                    l_c[i] = min(C(h_hat[i - 1] / (2c[i])), C(h_hat[i] / (-2c[i])))
                end

                return 1 .- l_c
            else
                return zeros(size(Rr, 1))
            end
        end

        if blend == :subcell
            ### Shock capturing optimization target
            b = sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations)) + sum(dot.(v, rhs_vol_low))        
            a = -dot.(vec(v' * Δr), (Rr * (rhs_vol_high - rhs_vol_low)))

            maximum_theta = preserve_positivity_forward()

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

            # determine blending factors
            θ = cache.knapsack_solver!(a, b, upper_bounds = maximum_theta * ones(length(a)))

            if nodewise_shock_capturing > 0
                @. a *= N^2 * tan(pi/2 * nodewise_shock_capturing)^2 * (θ - 1)^2 + 1

                θ = cache.knapsack_solver!(a, b)
            end

            view(du, :, e) .+= rd.M \ (Δr * (Diagonal(θ) * (Rr * rhs_vol_high) + 
                                            Diagonal(1 .- θ) * (Rr * rhs_vol_low)))
        elseif blend == :subcell_reversed
            # Get limiting coeffs
            l_c_1 = preserve_positivity_backward_nodewise(1)
            l_c_4 = preserve_positivity_backward_nodewise(4)
            # Get elementwise maxima
            l_c = max.(l_c_1, l_c_4)

            # Elementwise limiting coeffs for 2D sims
            if initial_condition == inital_condition_sedov_blastwave
                l_c .= maximum(l_c)
            end

            a = dot.(vec(v' * Δr), (Rr * (rhs_vol_low - rhs_vol_high)))
            b = -sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations)) - sum(dot.(v, rhs_vol_high)) - a'l_c

            # Call the Knapsack Solver
            θ = cache.knapsack_solver!(a, b, upper_bounds = (1 .- l_c))
            
            if nodewise_shock_capturing > 0
                @. θ *= N^(2nodewise_shock_capturing) * θ^2 / ((1 - θ)^2 + 1e-10) + 1
                θ = clamp.(θ, 0, 1 .- l_c)
            end
            
            view(du, :, e) .+= rd.M \ (rhs_vol_high + Δr * Diagonal(θ + l_c) * Rr * (rhs_vol_low - rhs_vol_high))
        elseif blend == :viscosity
            l = sum(dot.(v, rhs_vol_visc))
            dΨ = sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations))
            Ec = sum(dot.(v, rhs_vol_high))
            λ = 0

            if (0 > dΨ + Ec)
                a = 2(dΨ + Ec)
                b = l
                λ = a * b / (b^2 + 1e-12)
            end

            r = (rhs_vol_high - λ / 2 * rhs_vol_visc)

            # if !(-sum(dot.(v, r)) <= dΨ + 300 * eps())
            #     println("VIOLATION: ", -sum(dot.(v, r)) - dΨ)
            # end

            view(du, :, e) .+= rd.M \ r
        elseif blend == :hyperviscosity
            w1 = 1
            w2 = 1
            l_visc = sum(dot.(v, rhs_vol_visc))
            l_yimin = sum(dot.(v, rhs_vol_low))
            dΨ = sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations))
            Ec = sum(dot.(v, rhs_vol_high))
            b = dΨ + Ec

            a = [w1 * l_visc, w2 * l_yimin]

            θ = [0, 0]

            if a'a > 1e-12 && (b < -100 * eps())
                θ = b * a / (a'a + 1e-12)
            end

            r = (rhs_vol_high - w1 * θ[1] * rhs_vol_visc - w2 * θ[2] * rhs_vol_low)

            view(du, :, e) .+= rd.M \ r
        elseif blend == :nodewiseviscosity
            a = dot.(vec(v' * Δr), (Rr * (rhs_vol_visc)))
            ac = clamp.(a, 0, Inf)

            dΨ = sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations))
            Ec = sum(dot.(v, rhs_vol_high))
            λ = 0 * a

            if (0 > dΨ + Ec)
                b = dΨ + Ec
                λ = b * a / (a'a)
            end

            r = (rhs_vol_high - Δr * Diagonal(λ) * Rr * rhs_vol_visc)

            view(du, :, e) .+= rd.M \ r
        end

    end

    if any(isnan.(sum(du)))
        println("NAN DETECTED")
    end

    @. du /= -md.J


    return du
end

rd = RefElemData(Quad(), SBP(), N)
md = MeshData(uniform_mesh(rd.element_type, K), rd; 
is_periodic=true)

Qr = rd.M * rd.Dr 
Qs = rd.M * rd.Ds 
Qr_skew, Qs_skew = map(A -> (A - A'), (Qr, Qs))

(Qr_sparse, Qs_sparse), E_sparse = sparse_low_order_SBP_operators(rd)
Qr_sparse_skew, Qs_sparse_skew = map(A -> (A - A'), (Qr_sparse, Qs_sparse))

(Δr, Δs), (Rr, Rs) = subcell_limiting_operators(rd)
# override subcell limiting operators for preserve_positivity
Δr = diagm(0 => -ones(size(Rr, 2)), 1=>ones(size(Rr, 2)))[1:end-1,:]
Rr = tril(ones(size(Δr')), -1)

current_timestep = dt # for positivity preservation