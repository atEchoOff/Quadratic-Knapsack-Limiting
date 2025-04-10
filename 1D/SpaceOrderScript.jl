include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

using LsqFit

MODULE = "density_wave.jl"

N = 3
adaptive = false
dt = 1e-4

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). Last, N = $N. ENSURE CORRECT BEFORE PROCEEDING")

# sleep(15)

Li = Float64[]
Ki = Float64[]

i = 8
while length(Li) == 0 || Li[end] > 1e-6
    global K, i
    K = floor(Int, i)
    println("Testing K = $K")
    i *= 2^.2

    push!(Ki, 2 / K)

    include(MODULE)

    nm = L2_error
    println("Obtained norm $nm")
    push!(Li, nm)
end


println(Polynomials.fit(log.(Ki), log.(Li), 1))
display(Li)