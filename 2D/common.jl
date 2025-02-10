function rhs!(du, u, cache, t)
    (; rd, md) = cache
    (; Qr_skew, Qs_skew) = cache
    (; Qr_sparse_skew, Qs_sparse_skew, Δr, Rr) = cache
    (; equations, volume_flux, surface_flux) = cache
    (; VDM_inv, shock_capturing) = cache

    # interface flux
    uf = rd.Vf * u
    uP = uf[md.mapP]
    interface_flux = @. surface_flux(uf, uP, 
                                     SVector(md.nx, md.ny), equations) * md.Jf
    du .= rd.LIFT * interface_flux

    viscosities = Vector{Float64}(undef, K^2)
   
    @inbounds Threads.@threads for e in axes(du, 2)
    # for e in axes(du, 2)
        θ = zeros(size(Rr, 1))
        rhs_vol_high = similar(du, size(du, 1))
        rhs_vol_low = similar(rhs_vol_high)
        rhs_vol_visc = similar(rhs_vol_high)
        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))
        fill!(rhs_vol_visc, zero(eltype(rhs_vol_visc)))

        # if blend == :viscosity
        #     # Add surface contribution to r, and zero out du
        #     rhs_vol_high .= view(du, :, e)
        #     fill!(view(du, :, e), zero(eltype(du)))
        # end

        u_element = view(u, :, e)
        v = cons2entropy.(u_element, equations)
        v_avg = cons2entropy(.5 * sum(rd.wq .* u_element), equations)

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
                        Fij = flux_lax_friedrichs(u[i, e], u[j, e], nij, equations)
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

        if blend == :subcell
            ### Shock capturing optimization target
            b = sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations)) + sum(dot.(v, rhs_vol_low))        
            a = -dot.(vec(v' * Δr), (Rr * (rhs_vol_high - rhs_vol_low)))

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
            fill!(θ, one(eltype(θ)))
            θ = cache.knapsack_solver!(vec(a), b)

            if nodewise_shock_capturing > 0
                @. a *= N^2 * tan(pi/2 * nodewise_shock_capturing)^2 * (θ - 1)^2 + 1

                θ = cache.knapsack_solver!(vec(a), b)
            end

            view(du, :, e) .+= rd.M \ (Δr * (Diagonal(θ) * (Rr * rhs_vol_high) + 
                                            Diagonal(1 .- θ) * (Rr * rhs_vol_low)))
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

            viscosities[e] = λ / 2

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

    if blend == :viscosity
        # Add the surface term back
        interface_flux = ((uf - uP) .* md.nxJ) * Diagonal(viscosities)
        du .+= rd.LIFT * interface_flux
    end

    @. du /= -md.J
    
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