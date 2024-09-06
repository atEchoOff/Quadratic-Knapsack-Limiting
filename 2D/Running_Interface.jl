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

N = 3
K = 64

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))

# volume_flux = flux_ranocha
volume_flux = flux_central
# volume_flux = flux_shima_etal # useful for non ec solvers, for KHI

blend = :subcell

shock_capturing = false

abstol = 1e-6
reltol = 1e-4

timestepper = Tsit5()
adaptive = true
dt = 1e-4
saveat = .05