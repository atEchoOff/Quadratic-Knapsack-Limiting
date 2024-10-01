include("../ObjectTransitioner.jl")
using Plots

L1 = find("ContinuousKnapsackSolver")
L2 = find("QuadraticKnapsackSolver")
La = find("QuadraticKnapsackSolverA")
ES = find("NonKnapsackSolver")

# t = 2 .^ ((0:35) * .2) * 5e-6

plot(t, L1, lw=2, c=:red, label="L1", xaxis=:log, yaxis=:log)
plot!(t, L2, lw=2, c=:blue, label="L2", xaxis=:log, yaxis=:log)
plot!(t, La, lw=2, c=:purple, label="La", xaxis=:log, yaxis=:log)
plot!(t, ES, lw=2, c=:orange, label="ES", xaxis=:log, yaxis=:log)
xlabel!("dt")
ylabel!("L2 Error")
plot!(legend=:bottomright)