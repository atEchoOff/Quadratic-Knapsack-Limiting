include("common.jl")
tspan = (0, 0.2)
psi(u, ::CompressibleEulerEquations1D) = u[2] # rho * v1

function domain_change(x)
    a = 0.0
    b = 1.0
    return (b - a)/2 * (x + 1) + a
end

# Reset mesh data for domain change
VX = domain_change.(VX)
md = MeshData((VX, ), EToV, rd)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

equations = CompressibleEulerEquations1D(1.4)

initial_condition = initial_condition_modified_sod

u0 = rd.Pq * initial_condition.(md.xq, 0.0, equations)
u = copy(u0)

# Build the cache for RHS!
cache = (; 
           rd, md, 
           B = Diagonal([-1; zeros(rd.N-1); 1]),
           high_order_operators=(; Q_skew = rd.M * rd.Dr - (rd.M * rd.Dr)'), 
           low_order_operators=(; Q_skew_low = Qr-Qr'), 
           fv_operators = (; Δ, R),
           bc = [u0[1, 1] u0[end, end]],
           VDM_inv = inv(rd.VDM),
           shock_capturing = shock_capturing,
           nodewise_shock_capturing = nodewise_shock_capturing,
           blend = blend,
        knapsack_solver = knapsack_solver,
        physics = (; equations, volume_flux = volume_flux, surface_flux = flux_lax_friedrichs),
        );
ode = ODEProblem(rhs!, u, tspan, cache)

sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, callback=AliveCallback(alive_interval=1000), adaptive=adaptive, saveat=saveat)

println("Completed run with N = $N, K = $K, knapsack_solver = $(typeof(knapsack_solver)), timestepper = $(typeof(timestepper)), abstol = $abstol, reltol = $reltol")

u = sol.u[end]
u_plot = rd.Vp * getindex.(u, 1)
plot(rd.Vp * md.x, u_plot, leg=false)

# @gif for u in sol.u
#     plot(x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0, 1))
# end