include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

MODULE = "Modified_Sod_Shocktube.jl"

@assert adaptive == false
@assert dt == 1e-6

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). Last, (N, M) = ($N, $K). ENSURE CORRECT BEFORE PROCEEDING")

# sleep(15)
include(MODULE)

rq, wq = gauss_quad(0, 0, N+2)
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)

good = copy(sol.u[end])
Li = Float64[]
t = Float64[]

for i in (2 .^ ((0:35) * .2) * 5e-6) # modified sod
# for i in (2 .^ ((0:35) * .2) * 1e-5) # burgers
# for i in (2 .^ ((0:35) * .03) * 1e-3) # density wave
    global dt
    dt = i
    push!(t, dt)

    include(MODULE)

    nm = sqrt(sum(wJq .* map(x -> sum(x.^2), Vq * (good - sol.u[end]))))
    println("Obtained norm $nm")
    push!(Li, nm)
end

t_nonnan = t[isnan.(Li) .== false]
Li_nonnan = Li[isnan.(Li) .== false]
println(Polynomials.fit(log.(t_nonnan), log.(Li_nonnan), 1))

# Save all methods to object transfer
save("dts", t)
save(case, Li)