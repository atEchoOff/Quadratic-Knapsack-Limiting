using LinearAlgebra, SparseArrays
using StartUpDG
using Trixi
using OrdinaryDiffEq
using Plots

include("../L1_knapsack.jl")
include("../L2_knapsack.jl")

rd = RefElemData(Quad(), SBP(), N)
md = MeshData(uniform_mesh(rd.element_type, K), rd; 
              is_periodic=true)

Qr = rd.M * rd.Dr 
Qs = rd.M * rd.Ds 
Qr_skew, Qs_skew = map(A -> (A - A'), (Qr, Qs))

(Qr_sparse, Qs_sparse), E_sparse = sparse_low_order_SBP_operators(rd)
Qr_sparse_skew, Qs_sparse_skew = map(A -> (A - A'), (Qr_sparse, Qs_sparse))

(Δr, Δs), (Rr, Rs) = subcell_limiting_operators(rd)

# physics
equations = CompressibleEulerEquations2D(1.4)
# volume_flux = flux_ranocha
# volume_flux = flux_central
volume_flux = flux_shima_etal
surface_flux = flux_lax_friedrichs

cache = (; rd, md, Qr_skew, Qs_skew, 
           Qr_sparse_skew, Qs_sparse_skew, Δr, Rr, 
           equations, volume_flux, surface_flux,
        knapsack_solver! = knapsack_solver, VDM_inv=inv(rd.VDM), shock_capturing=shock_capturing)

# du/dt + dfx(u)/dx + dfy(u)/dy = 0
function initial_condition(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    # rho = 1 + 0.5 * (abs(x) < 0.5) * (abs(y) < 0.5)
    # v1 = 0.
    # v2 = 0.
    # p = rho^equations.gamma
    # return prim2cons(SVector(rho, v1, v2, p), equations)

    # KHI
    B = tanh(15 * y + 7.5) - tanh(15 * y - 7.5)
    rho = 0.5 + 0.75 * B
    u = 0.5 * (B - 1)
    v = 0.1 * sin(2 * pi * x)
    p = 1
    return prim2cons(SVector(rho, u, v, p), equations)
end

u0 = initial_condition.(SVector.(md.x, md.y), 0.0, equations)
du = similar(u0)

function psi(u, normal, equations::CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, _ = u
    return dot(SVector(rho_v1, rho_v2), normal)
end

function rhs!(du, u, cache, t)
    (; rd, md) = cache
    (; Qr_skew, Qs_skew) = cache
    (; Qr_sparse_skew, Qs_sparse_skew, Δr, Rr) = cache
    (; equations, volume_flux, surface_flux) = cache
    (; VDM_inv, shock_capturing) = cache

    # interface flux
    uf = rd.Vf * u
    interface_flux = @. surface_flux(uf, uf[md.mapP], 
                                     SVector(md.nx, md.ny), equations) * md.Jf
    du .= rd.LIFT * interface_flux
   
    θ = zeros(size(Rr, 1))
    rhs_vol_high = similar(du, size(du, 1))
    rhs_vol_low = similar(rhs_vol_high)
    for e in axes(du, 2)
        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))

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

                # low order; create a unit nij vector
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
        end

        # perform knapsack limiting - this can be made more efficient
        u_element = view(u, :, e)
        v = cons2entropy.(u_element, equations)
        b = sum(rd.wf .* psi.(uf[:, e], SVector.(md.nxJ[:,e], md.nyJ[:,e]), equations)) + sum(dot.(v, rhs_vol_low))        
        a = -dot.(vec(v' * Δr), (Rr * (rhs_vol_high - rhs_vol_low)))
        if !(b > -100 * eps()) # assert rhs is positive
            @warn "b is not positive"
            @show b
        end

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

        # determine blending factors
        fill!(θ, one(eltype(θ)))
        θ = cache.knapsack_solver!(vec(a), b)
        view(du, :, e) .+= rd.M \ (Δr * (Diagonal(θ) * (Rr * rhs_vol_high) + 
                                         Diagonal(1 .- θ) * (Rr * rhs_vol_low)))

    end

    @. du /= -md.J
    
end

# # check entropy conservation
# rhs!(du, u0, cache, 0.0)
# @show sum(dot.(cons2entropy.(u0, equations), rd.M * (du .* md.J)))

tspan = (0.0, 3.5)
ode = ODEProblem(rhs!, u0, tspan, cache)

sol = solve(ode, timestepper, 
            abstol=abstol, reltol=reltol,
            saveat=LinRange(tspan..., 25),
            callback=AliveCallback(alive_interval=10), adaptive=adaptive, dt=dt)

# interpolate solution to equispaced points for plotting
rp, sp = equi_nodes(rd.element_type, 5)
Vp = vandermonde(rd.element_type, N, rp, sp) / rd.VDM
u = sol.u[end]
scatter(vec(Vp * md.x), vec(Vp * md.y), 
            zcolor=vec(Vp * getindex.(u, 1)), 
            msw=0, ms=2, legend=false, ratio=1)

@gif for u in sol.u
    scatter(vec(Vp * md.x), vec(Vp * md.y), 
            zcolor=vec(Vp * getindex.(u, 1)), 
            msw=0, ms=2, legend=false, ratio=1)
end fps=10