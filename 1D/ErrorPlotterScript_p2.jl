include("../ObjectTransitioner.jl")
using Plots
using LaTeXStrings
using RollingFunctions

println(@__DIR__)

N_to_plot = 3
K_to_plot = 64
TS_to_plot = "Tsit5"

L1 = find("N$(N_to_plot)K$(K_to_plot)TS$(TS_to_plot)CL1")
L2 = find("N$(N_to_plot)K$(K_to_plot)TS$(TS_to_plot)CL2")
EC = find("N$(N_to_plot)K$(K_to_plot)TS$(TS_to_plot)CEC")
DG = find("N$(N_to_plot)K$(K_to_plot)TS$(TS_to_plot)CDG")
# times = find("times")
dt = 1e-5
T = .2
times = dt:dt:T

times = rollmean(times, 1024)
L1 = rollmean(L1, 1024)
L2 = rollmean(L2, 1024)
EC = rollmean(EC, 1024)
DG = rollmean(DG, 1024)

# times = times[2:end]

# times should already be set
# plot()
plot(times, L1, lw=2, c=:red, label="Linear Knapsack", yaxis=:log)
plot!(times, L2, lw=2, c=:blue, label="Quadratic Knapsack", yaxis=:log)
plot!(times, EC, lw=2, c=:orange, label="Entropy Conservative", yaxis=:log)
plot!(times, DG, lw=2, c=:green, label="DG", yaxis=:log)
xlabel!(L"$t$")
ylabel!("Total Entropy")
plot!(legend=:topleft)
plot!(legendfontsize=10)
plot!(xtickfontsize=10)
plot!(ytickfontsize=10)
plot!(labelfontsize=12)

# plot!(ylim=(10^-.65, 10^-.55))
# plot!(xlim=(.15, .2))

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