include("../ObjectTransitioner.jl")
using Plots

L1 = find("ContinuousKnapsackSolver")
L2 = find("QuadraticKnapsackSolver")
La = find("QuadraticKnapsackSolverA")
ES = find("NonKnapsackSolver")

# times should already be set

plot(times, L1, lw=2, c=:red, label="L1", yaxis=:log)
plot!(times, L2, lw=2, c=:blue, label="L2", yaxis=:log)
plot!(times, La, lw=2, c=:purple, label="La", yaxis=:log)
plot!(times, ES, lw=2, c=:orange, label="ES", yaxis=:log)
xlabel!("t")
ylabel!("L2 Error")
plot!(legend=:bottomright)