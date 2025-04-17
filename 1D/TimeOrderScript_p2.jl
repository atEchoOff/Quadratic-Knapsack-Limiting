include("../ObjectTransitioner.jl")
using Plots
using LaTeXStrings

# L1 = find("L1")
# L2 = find("L2")
# EC = find("EC")
# t = find("dts")

# t = 2 .^ ((0:35) * .2) * 5e-6

plot(t, L1, lw=2, c=:red, label="LK", xaxis=:log, yaxis=:log)
plot!(t, L2, lw=2, c=:blue, label="QK", xaxis=:log, yaxis=:log)
plot!(t, ES, lw=2, c=:orange, label="ESFD", xaxis=:log, yaxis=:log)
xlabel!(L"$\Delta t$")
ylabel!(L"$L^2$ Error")
niceplot!()
plot!(legend=:bottomright)