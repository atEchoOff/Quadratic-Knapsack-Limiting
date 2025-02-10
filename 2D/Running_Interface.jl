using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots

include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_weighted_a.jl")
include("../Non_knapsack.jl")
include("initial_conditions.jl")
include("../run_saver.jl")

dimstring = "2D"
use_run_saver = false

N = 4
K = 100

total_error_estimates = []

# if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
#     knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))
# end
knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))

# volume_flux = flux_ranocha
volume_flux = flux_central
# volume_flux = flux_shima_etal # useful for non ec solvers, for KHI

blend = :subcell

shock_capturing = 0
nodewise_shock_capturing = 0

abstol = 1e-6
reltol = 1e-4

timestepper = SSPRK43()
adaptive = true
dt = 1e-4
saveat = 1
