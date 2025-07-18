include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

MODULE = "Modified_Sod_Shocktube.jl"

# @assert adaptive == false
# @assert dt == 1e-5

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). Last, (N, M) = ($N, $K). ENSURE CORRECT BEFORE PROCEEDING")

# sleep(15)

include(MODULE)
times = copy(sol.t)

errors = Float64[]

rq, wq = gauss_quad(0, 0, N+2)
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x


for i in 1:length(times)
    # push!(errors, sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(xq, sol.t[i], equations) - Vq * sol.u[i])))) # for error
    push!(errors, md.J[1, 1] * sum(rd.wq .* (entropy.(sol.u[i], equations)))) # for total entropy
end

times = times[2:end]
errors = errors[2:end]

save(case, errors)
