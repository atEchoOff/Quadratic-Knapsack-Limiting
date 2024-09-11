include("../ObjectTransitioner.jl")
using Plots

L1 = find("ContinuousKnapsackSolver")
L2 = find("QuadraticKnapsackSolver")
La = find("QuadraticKnapsackSolverA")
ES = find("NonKnapsackSolver")

# times = times[2:end]

# times should already be set

plot(times, L1, lw=2, c=:red, label="L1", yaxis=:log)
plot!(times, L2, lw=2, c=:blue, label="L2", yaxis=:log)
plot!(times, La, lw=2, c=:purple, label="La", yaxis=:log)
plot!(times, ES, lw=2, c=:orange, label="ES", yaxis=:log)
xlabel!("t")
ylabel!("L2 Error")
plot!(legend=:topleft)

# savefig("temp.png")

# plot(times[times .>= .15], L1[times .>= .15], lw=2, c=:red, label="L1", yaxis=:log)
# plot!(times[times .>= .15], L2[times .>= .15], lw=2, c=:blue, label="L2", yaxis=:log)
# # plot!(times[times .>= .15], La[times .>= .15], lw=2, c=:purple, label="La", yaxis=:log)
# plot!(times[times .>= .15], ES[times .>= .15], lw=2, c=:orange, label="ES", yaxis=:log)
# xlabel!("t")
# ylabel!("L2 Error")
# plot!(legend=:topleft)

# savefig("temp_zoomed_at_end.png")

# plot(rollmean(times, 1024), rollmean(L1, 1024), lw=2, c=:red, label="L1", yaxis=:log)
# plot!(rollmean(times, 1024), rollmean(L2, 1024), lw=2, c=:blue, label="L2", yaxis=:log)
# # plot!(rollmean(times, 1024), rollmean(La, 1024), lw=2, c=:purple, label="La", yaxis=:log)
# plot!(rollmean(times, 1024), rollmean(ES, 1024), lw=2, c=:orange, label="ES", yaxis=:log)
# xlabel!("t")
# ylabel!("L2 Error")
# plot!(legend=:topleft)

# savefig("temp_rolled.png")