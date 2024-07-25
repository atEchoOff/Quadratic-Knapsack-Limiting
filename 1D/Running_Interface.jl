using Pkg

Pkg.activate(".")

using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots
using Polynomials

# include("../Space_Error_Calculator.jl")
include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_weighted_a.jl")
include("../Non_knapsack.jl")
include("../Smooth_Knapsack.jl")

N = 3
K = 4

include("initial_conditions.jl")

# knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
knapsack_solver = ContinuousKnapsackSolver((N + 2))
# knapsack_solver = SmoothKnapsackSolver{Float64}()

volume_flux = flux_central
# volume_flux = flux_ranocha
# volume_flux = flux_ec
blend = :subcell

shock_capturing = false
nodewise_shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = Tsit5()
adaptive = true
dt = 1e-6

total_error_estimates = Float64[]

# include("Modified_Sod_Shocktube.jl")
# include("Modified_Sod_Shocktube_Dynamic_Theta.jl")
# include("density_wave.jl")
include("burgers.jl")
# include("advection.jl")


### TIME CONVERGENCE
# good = copy(sol.u[end])
# Li = Float64[]
# t = Float64[]

# for i in 2 .^ ((0:35) * .2)
#     global dt
#     dt = i * 1e-5
#     push!(t, dt)

#     include("burgers.jl")

#     nm = sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), good - sol.u[end])))
#     println("Obtained norm $nm")
#     push!(Li, nm)
# end

# println(fit(log.(Ki), log.(Li), 1))

### SPACE CONVERGENCE
Li = Float64[]
Ki = Float64[]

for i in 2 .^ ((8:20) * .5)
    global K
    K = floor(Int, i)
    push!(Ki, 1 / K)

    include("burgers.jl")

    nm = L2_error
    println("Obtained norm $nm")
    push!(Li, nm)
end

println(fit(log.(Ki), log.(Li), 1))