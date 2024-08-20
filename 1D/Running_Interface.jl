using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots
using Polynomials

include("../CustomIntegrator.jl")
include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_weighted_a.jl")
include("../Non_knapsack.jl")
include("../Smooth_Knapsack.jl")
include("initial_conditions.jl")

N = 4
K = 256

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 2))
# knapsack_solver = SmoothKnapsackSolver{Float64}()

total_error_estimates = []

volume_flux = flux_central
# volume_flux = flux_ranocha
# volume_flux = flux_ec
blend = :subcell

shock_capturing = false
nodewise_shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = SSPRK43()
adaptive = true
dt = 1e-4
saveat = .05
