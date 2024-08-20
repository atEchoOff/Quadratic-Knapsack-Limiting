include("common.jl")
tspan = (0.0, 8.0)
psi(u, ::InviscidBurgersEquation1D) = (1/6 * u .^ 3)[1]

md = make_periodic(md)

equations = InviscidBurgersEquation1D()

initial_condition = initial_condition_burgers

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
@time sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, saveat=saveat, callback=AliveCallback(alive_interval=1000), adaptive=adaptive)

# @gif for u in sol.u
#     plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0, 2))
# end
u = sol.u[end]
# plot!(rd.Vp * md.x, getindex.(initial_condition.(rd.Vp * md.x, tspan[2], equations), 1), leg=false)

rq, wq = gauss_quad(0, 0, N+2)
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x

function true_sol(x)
    # return SVector{1}(.25 * x + .5)
    return SVector{1}(.2 * x)
end

L1_error = sum(wJq .* map(x -> sum(abs.(x)), true_sol.(xq) - Vq * sol.u[end]))
L2_error = sqrt(sum(wJq .* map(x -> sum(x.^2), true_sol.(xq) - Vq * sol.u[end])))

println("N = $N, K1D = $(md.num_elements), L1_error = $L1_error, L2_error = $L2_error")


plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false)

@gif for u in sol.u
    plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0., 1.2))
end