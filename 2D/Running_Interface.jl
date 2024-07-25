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

N = 3
K = 64

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))
blend = :subcell

total_error_estimates = Float64[]

shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = Tsit5()
adaptive = true
dt = 1e-4

# include("KHI.jl")
include("density_wave.jl")
# include("Riemann.jl")