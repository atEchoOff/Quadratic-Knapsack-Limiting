# Instructions:
# Make sure there is no saveat in MST
# First, include custom CustomIntegrator.jl file.
# Then, run with adaptive = true to get the tolerances saved in total_error_tolerances
# Save solutions = copy(sol.u)
# Save tolerances = copy(total_error_estimates)
# Save times = copy(sol.t)
# Now, edit Modified_Sod_Shocktube.jl to have saveat = times
# Set adaptive = false with a small dt, and run MST again
# Run this to get errors, and plot errors against tolerances
MODULE = "Modified_Sod_Shocktube.jl"

# Less refined solution
include(MODULE)
solutions = copy(sol.u)
times = copy(sol.t)

@assert all(isnan.(getindex.(sol.u[end], 1)) .== false)
# saveat = times

# Refined solution
adaptive = false
dt = 1e-6
include(MODULE)

errors = []
for i in 1:length(times)
    push!(errors, sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), solutions[i] - sol.u[i]))))
end

times = times[2:end]
errors = errors[2:end]
