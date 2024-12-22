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

N = 10
K = 64

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
# volume_flux = flux_shima_etal
# volume_flux = flux_lax_friedrichs
# volume_flux = flux_ranocha
# volume_flux = flux_ec
blend = :hyperviscosity

shock_capturing = 0
nodewise_shock_capturing = 0

abstol = 1e-8
reltol = 1e-6

timestepper = Tsit5()
adaptive = true
dt = 1e-4
saveat = []
