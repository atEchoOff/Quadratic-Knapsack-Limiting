using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots

include("../CustomIntegrator.jl")
include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_minimizer.jl")
include("../L2_knapsack_weighted_a.jl")
include("../Non_knapsack.jl")
include("initial_conditions.jl")
include("../run_saver.jl")

dimstring = "2D"
use_run_saver = false

N = 3
K = 32

total_error_estimates = []

knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackMinimizer{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))

# volume_flux = flux_ranocha
volume_flux = flux_central
# volume_flux = flux_shima_etal # useful for non ec solvers, for KHI

blend = :subcell

shock_capturing = 0
nodewise_shock_capturing = 0

abstol = 1e-8
reltol = 1e-6

timestepper = SSPRK43()
adaptive = false
dt = 2e-5
saveat = .1

preserve_positivity = .8

# for automation purposes
if @isdefined case
    if case == "L1"
        global knapsack_solver = ContinuousKnapsackSolver((N + 1) * (N + 2))
        global volume_flux = flux_central
    elseif case == "L2"
        global knapsack_solver = QuadraticKnapsackSolver{Float64}()
        global volume_flux = flux_central
    elseif case == "EC"
        global knapsack_solver = NonKnapsackSolver{Float64}()
        global volume_flux = flux_ranocha
    elseif case == "DG"
        global knapsack_solver = NonKnapsackSolver{Float64}()
        global volume_flux = flux_central
    end
end