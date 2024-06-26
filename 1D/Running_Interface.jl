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

N = 7
K = 128

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver(N + 2)
blend = :subcell

shock_capturing = false

abstol = 1e-8
reltol = 1e-6

timestepper = RDPK3SpFSAL35()
adaptive = true
dt = 5e-5

include("Modified_Sod_Shocktube.jl")
# include("limited_entropy_stable_tests_jesse.jl")
# include("Global_Entropy_Stable_Tests.jl")
# include("Global_Modified_Sod_Shocktube.jl")