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

N = 4
K = 128

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))
blend = :subcell

shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = SSPRK43()
adaptive = true
dt = 1e-4

# include("KHI.jl")
# include("Advection.jl")
include("Riemann.jl")