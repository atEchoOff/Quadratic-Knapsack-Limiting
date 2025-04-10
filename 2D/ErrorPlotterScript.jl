include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

MODULE = "density_wave.jl"

@assert adaptive == false

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). Last, (N, M) = ($N, $K). ENSURE CORRECT BEFORE PROCEEDING")

include(MODULE)
times = copy(sol.t)
errors = Float64[]

rq, wq = gauss_quad(0, 0, N+2)
sq = copy(rq)
Vq = vandermonde(rd.element_type, rd.N, rq, sq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x
yq = Vq * md.y

for i in 1:length(times)
    push!(errors, sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(SVector.(xq, yq), times[i], equations) - Vq * sol.u[i]))))
    # push!(errors, sum(md.wJq .* (Vq * entropy.(sol.u[i], equations))))
end

times = times[2:end]
errors = errors[2:end]

save(case, errors)