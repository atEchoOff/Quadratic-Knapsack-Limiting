include("common.jl")
tspan = (0.0, 10.)
equations = CompressibleEulerEquations2D(1.4)

cache = (; rd, md, Qr_skew, Qs_skew, 
           Qr_sparse_skew, Qs_sparse_skew, Δr, Rr, 
           equations, volume_flux, surface_flux = flux_lax_friedrichs,
        knapsack_solver! = knapsack_solver, VDM_inv=inv(rd.VDM), shock_capturing=shock_capturing)

initial_condition = initial_condition_khi

u0 = initial_condition.(SVector.(md.x, md.y), 0.0, equations)
du = similar(u0)

function psi(u, normal, equations::CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, _ = u
    return dot(SVector(rho_v1, rho_v2), normal)
end

# # check entropy conservation
# rhs!(du, u0, cache, 0.0)
# @show sum(dot.(cons2entropy.(u0, equations), rd.M * (du .* md.J)))
ode = ODEProblem(rhs!, u0, tspan, cache)

sol = load_save_if_exists()
if sol == -1
    sol = solve(ode, timestepper, 
            abstol=abstol, reltol=reltol,
            saveat=saveat,
            callback=AliveCallback(alive_interval=10), adaptive=adaptive, dt=dt)
    save_run()
end

# interpolate solution to equispaced points for plotting
rp, sp = equi_nodes(rd.element_type, 5)
Vp = vandermonde(rd.element_type, N, rp, sp) / rd.VDM
u = sol.u[end]
scatter(vec(Vp * md.x), vec(Vp * md.y), 
            zcolor=vec(Vp * getindex.(u, 1)), 
            msw=0, ms=2, legend=false, ratio=1)

# @gif for u in sol.u
#     scatter(vec(Vp * md.x), vec(Vp * md.y), 
#             zcolor=vec(Vp * getindex.(u, 1)), 
#             msw=0, ms=2, legend=false, ratio=1)
# end fps=10
