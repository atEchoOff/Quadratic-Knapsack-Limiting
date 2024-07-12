using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots
using CSV
using Tables

include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_weighted_a.jl")
include("../Non_knapsack.jl")
include("../Smooth_Knapsack.jl")

N = 7
K = 32

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver(N + 2)
# knapsack_solver = SmoothKnapsackSolver{Float64}()

volume_flux = flux_central
# volume_flux = flux_ranocha
# volume_flux = flux_ec
blend = :subcell

shock_capturing = false
nodewise_shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = RK4()
adaptive = true
dt = 1e-5

total_error_estimates = Float64[]

# include("Stationary_Contact.jl")
include("Modified_Sod_Shocktube.jl")
# include("Modified_Sod_No_Low_Order.jl")
# include("limited_entropy_stable_tests_jesse.jl")
# include("burgers.jl")
# include("advection.jl")
# include("Global_Entropy_Stable_Tests.jl")
# include("Global_Modified_Sod_Shocktube.jl")

# good = copy(sol.u[end])
# Li = Float64[]
# t = Float64[]

# for i in 2 .^ ((0:35) * .2)
#     global dt
#     dt = i * 1e-4
#     push!(t, dt)

#     include("burgers.jl")

#     nm = sqrt(sum(diag(rd.M) .* map(x -> sum(x.^2), good - sol.u[end])))
#     println("Obtained norm $nm")
#     push!(Li, nm)
# end

# p0 = [.5, .5, .5]
# model(t, p) = p[1] .+ p[2] * t.^(p[3])
# println(fit(log.(t), log.(Li), 1))
