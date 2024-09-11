include("../Jacobiator.jl")
include("../ObjectTransitioner.jl")
include("Running_Interface.jl")

MODULE = "burgers.jl"

@assert adaptive == false
@assert dt == 1e-5

println("Running $MODULE, with knapsack_solver = $(nameof(typeof(knapsack_solver))), volume_flux = $(nameof(typeof(volume_flux))) timestepper = $(nameof(typeof(timestepper))). ENSURE CORRECT BEFORE PROCEEDING")

# sleep(15)

# Override solve to instead compute maximal real eigenvalue
solve_old = solve

max_eigenvalues = Float64[]

function solve_new(ode, timestepper; args...)
    du = similar(u0)
    max_eigenvalue = maximum(real.(eigvals(get_jacobian_rhs(du, u0))))

    push!(max_eigenvalues, max_eigenvalue)
    error("Thrown only to stop testing in modules")
end

solve = solve_new

N = 3
K = 64

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

try
    include(MODULE)
catch
end

N = 7
K = 32

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

try
    include(MODULE)
catch
end

N = 15
K = 16

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

try
    include(MODULE)
catch
end

solve = solve_old

max_eigenvalues