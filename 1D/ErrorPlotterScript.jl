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
for i in 1:length(times)
    push!(errors, sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), initial_condition.(md.xq, times[i], equations) - sol.u[i]))))
    # push!(errors, sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), entropy.(sol.u[i], equations)))))
end

times = times[2:end]
errors = errors[2:end]

save(nameof(typeof(knapsack_solver)), errors)