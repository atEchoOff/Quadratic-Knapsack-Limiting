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

N = 10
K = 16

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver(N + 2)
blend = :subcell

shock_capturing = false
nodewise_shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = Midpoint()
adaptive = false
dt = 5e-6

# include("Modified_Sod_Shocktube.jl")
# include("limited_entropy_stable_tests_jesse.jl")
# include("Global_Entropy_Stable_Tests.jl")
# include("Global_Modified_Sod_Shocktube.jl")

for i in 1:20
    global dt
    dt = i * 1e-5

    include("Modified_Sod_Shocktube.jl")

    nm = norm(sol.u[end] - good_L2)
    println("Obtained norm $nm")
    push!(L2, nm)
end