include("common.jl")
tspan = (0, 10.)
psi(u, ::LinearScalarAdvectionEquation1D) = .5 * u^2

md = make_periodic(md)

equations = LinearScalarAdvectionEquation1D(1.0)

initial_condition = initial_condition_linear_advection

u0 = rd.Pq * initial_condition.(md.xq, 0.0, equations)
u = copy(u0)

cache = (; 
           rd, md, 
           B = Diagonal([-1; zeros(rd.N-1); 1]),
           high_order_operators=(; Q_skew = rd.M * rd.Dr - (rd.M * rd.Dr)'), 
           low_order_operators=(; Q_skew_low = Qr-Qr'), 
           fv_operators = (; Δ, R), 
           blend = blend,
           knapsack_solver = knapsack_solver,
           bc = nothing,
        physics = (; equations, volume_flux = volume_flux, surface_flux = flux_lax_friedrichs),
        shock_capturing = shock_capturing,
        nodewise_shock_capturing = nodewise_shock_capturing,
        VDM_inv = inv(rd.VDM)
        )

# dt = 0.25 * estimate_h(rd, md) * (1 / (2 * rd.N + 1))
ode = ODEProblem(rhs!, u, tspan, cache)

sol = load_save_if_exists()
if sol == -1
    sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, saveat=saveat, callback=AliveCallback(alive_interval=1000), adaptive=adaptive)
    save_run()
end
@gif for u in sol.u
    plot(md.x, getindex.(u, 1), leg=false, ylims=(0, 2))
end
# u = sol.u[end]
# # plot!(rd.Vp * md.x, getindex.(initial_condition.(rd.Vp * md.x, tspan[2], equations), 1), leg=false)

# rq, wq = gauss_quad(0, 0, N+2)
# Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
# wJq = Diagonal(wq) * (Vq * md.J)
# xq = Vq * md.x

# L1_error = sum(wJq .* map(x -> sum(abs.(x)), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end]))
# L2_error = sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end])))

# println("N = $N, K1D = $(md.num_elements), L1_error = $L1_error, L2_error = $L2_error")


# plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false)
