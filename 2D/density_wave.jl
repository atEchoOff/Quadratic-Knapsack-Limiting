include("common.jl")
tspan = (0.0, 4.0)
equations = CompressibleEulerEquations2D(1.4)

cache = (; rd, md, Qr_skew, Qs_skew, 
           Qr_sparse_skew, Qs_sparse_skew, Δr, Rr, 
           equations, volume_flux, surface_flux = flux_lax_friedrichs,
        knapsack_solver! = knapsack_solver, VDM_inv=inv(rd.VDM), shock_capturing=shock_capturing)

initial_condition = initial_condition_density_wave

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

sol = solve(ode, timestepper, 
            abstol=abstol, reltol=reltol,
            saveat=saveat,
            callback=AliveCallback(alive_interval=10), adaptive=adaptive, dt=dt)

rq, wq = gauss_quad(0, 0, N+2)
sq = copy(rq)
Vq = vandermonde(rd.element_type, rd.N, rq, sq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x
yq = Vq * md.y

L1_error = sum(wJq .* map(x -> sum(abs.(x)), initial_condition.(SVector.(xq, yq), sol.t[end], equations) - Vq * sol.u[end]))
L2_error = sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(SVector.(xq, yq), sol.t[end], equations) - Vq * sol.u[end])))

println("N = $N, K1D = $(md.num_elements), L1_error = $L1_error, L2_error = $L2_error")