MODULE = "Modified_Sod_Shocktube.jl"

include("Running_Interface.jl")
include(MODULE)

good = copy(sol.u[end])
Li = Float64[]
t = Float64[]

for i in 2 .^ ((0:35) * .2)
    global dt
    dt = i * 1e-5
    push!(t, dt)

    include(MODULE)

    nm = sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), good - sol.u[end])))
    println("Obtained norm $nm")
    push!(Li, nm)
end

println(fit(log.(t), log.(Li), 1))