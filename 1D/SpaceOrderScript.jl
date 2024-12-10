include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

MODULE = "density_wave.jl"

adaptive = false
dt = 1e-4

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). Last, N = $N. ENSURE CORRECT BEFORE PROCEEDING")

# sleep(15)

Li = Float64[]
Ki = Float64[]

for i in floor.(Int64, 2 .^ ((40:52) * .1))
    global K
    K = i
    push!(Ki, 1 / K)

    include(MODULE)

    nm = sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), initial_condition.(md.xq, sol.t[end], equations) - sol.u[end])))
    println("Obtained norm $nm")
    push!(Li, nm)
end

println(Polynomials.fit(log.(Ki), log.(Li), 1))