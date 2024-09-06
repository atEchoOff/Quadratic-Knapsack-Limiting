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

N = 7
K = 32

# if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
#     knapsack_solver = ContinuousKnapsackSolver((N + 2))
# end
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

timestepper = RK4()
adaptive = false
dt = 1e-5
saveat = .05
