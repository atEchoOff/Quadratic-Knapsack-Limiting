include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

MODULE = "burgers.jl"

@assert adaptive == false
@assert dt == 5e-6

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). Last, (N, M) = ($N, $K). ENSURE CORRECT BEFORE PROCEEDING")

sleep(15)
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

t_nonnan = t[isnan.(Li) .== false]
Li_nonnan = Li[isnan.(Li) .== false]
println(Polynomials.fit(log.(t_nonnan), log.(Li_nonnan), 1))

# Save all methods to object transfer
save(nameof(typeof(knapsack_solver)), Li)