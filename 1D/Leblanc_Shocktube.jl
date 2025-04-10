include("common.jl")
tspan = (0, .0001)
psi(u, ::CompressibleEulerEquations1D) = u[2] # rho * v1

function domain_change(x)
    a = -10
    b = 10
    return (b - a)/2 * (x + 1) + a
end

# Reset mesh data for domain change
VX = domain_change.(VX)
md = MeshData((VX, ), EToV, rd)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

equations = CompressibleEulerEquations1D(1.4)

initial_condition = initial_condition_leblanc_shocktube

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
        physics = (; equations, volume_flux = volume_flux, surface_flux = flux_hllc),
        );
ode = ODEProblem(rhs!, u, tspan, cache)

sol = load_save_if_exists()
if sol == -1
    sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, callback=AliveCallback(alive_interval=100), adaptive=adaptive, saveat=saveat)
    save_run()
end

println("Completed run with N = $N, K = $K, knapsack_solver = $(typeof(knapsack_solver)), timestepper = $(typeof(timestepper)), abstol = $abstol, reltol = $reltol")

u = sol.u[end]
u .= cons2prim.(u, equations)
u_plot = rd.Vp * getindex.(u, 1)
plot(vec(md.x), vec(getindex.(u, 1)), lw=2, leg=false, yaxis=:log)

# @gif for u in sol.u
#     plot(x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0, 1))
# end