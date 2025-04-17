include("../ObjectTransitioner.jl")
using Plots
using LaTeXStrings
using RollingFunctions
using Measures

println(@__DIR__)

L1 = find("L1")
L2 = find("L2")
EC = find("ES")
DG = find("DG")
# times = find("times")
dt = 1e-2
T = 4
times = dt:dt:T

# times = rollmean(times, 1024)
# L1 = rollmean(L1, 1024)
# L2 = rollmean(L2, 1024)
# EC = rollmean(EC, 1024)
# DG = rollmean(DG, 1024)

# times = times[2:end]

# times should already be set
plot()
plot(times, L1, lw=2, c=:red, label="LK", yaxis=:log)
plot!(times, L2, lw=2, c=:blue, label="QK", yaxis=:log)
plot!(times, EC, lw=2, c=:orange, label="ESFD", yaxis=:log)
plot!(times, DG, lw=2, c=:green, label="DGSEM", yaxis=:log)
# plot(times, L1, lw=2, c=:red, label="LK")
# plot!(times, L2, lw=2, c=:blue, label="QK")
# plot!(times, EC, lw=2, c=:orange, label="EC")
# plot!(times, DG, lw=2, c=:green, label="DG")
xlabel!(L"$t$")
ylabel!(L"$L^2$ Error")
plot!(legend=:right)
plot!(legendfontsize=12)
plot!(xtickfontsize=12)
plot!(ytickfontsize=12)
plot!(labelfontsize=14)
# plot!(rightmargin = 4mm)

L1 = L1[times .>= .1]
L2 = L2[times .>= .1]
EC = EC[times .>= .1]
DG = DG[times .>= .1]



plot!(ylim=(minimum(union(L1, L2, EC, DG)), maximum(union(L1, L2, EC, DG))))

savefig("density_wave_plots.png")
plot!()
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