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
include("../L2_knapsack_minimizer.jl")
include("../L2_knapsack_weighted_a.jl")
include("../Non_knapsack.jl")
include("../Smooth_Knapsack.jl")
include("initial_conditions.jl")
include("../run_saver.jl")

dimstring = "1D"
use_run_saver = false

N = 3
K = 100

knapsack_stats = Nothing # comment this out if you dont want to save blending coeffs and iters over time
# knapsack_solver = QuadraticKnapsackSolver{Float64}()
# knapsack_solver = NonKnapsackSolver{Float64}()
# knapsack_solver = QuadraticKnapsackSolverA{Float64}()
# knapsack_solver = ContinuousKnapsackSolver((N + 2))
# knapsack_solver = SmoothKnapsackSolver{Float64}()
knapsack_solver = QuadraticKnapsackMinimizer{Float64}()

total_error_estimates = []

volume_flux = flux_central
# volume_flux = flux_shima_etal
# volume_flux = flux_lax_friedrichs
# volume_flux = flux_ranocha
# volume_flux = flux_ec
blend = :subcell_reversed

shock_capturing = 0
nodewise_shock_capturing = 0

abstol = 1e-6
reltol = 1e-4

timestepper = RK4()
adaptive = true
dt = 1e-4
saveat = 1e-2

preserve_positivity = -1

# for automation purposes
if @isdefined case
    if case == "L1"
        global knapsack_solver = ContinuousKnapsackSolver((N + 2))
        global volume_flux = flux_central
        global blend = :subcell
    elseif case == "L2"
        global knapsack_solver = QuadraticKnapsackMinimizer{Float64}()
        global volume_flux = flux_central
        global blend = :subcell_reversed
    elseif case == "ES"
        global knapsack_solver = NonKnapsackSolver{Float64}()
        global volume_flux = flux_ranocha
        global blend = :subcell
    elseif case == "DG"
        global knapsack_solver = NonKnapsackSolver{Float64}()
        global volume_flux = flux_central
        global blend = :subcell
    end
end