MODULE = "burgers.jl"

println("Running $MODULE, with knapsack_solver = $(typeof(knapsack_solver)), volume_flux = $(typeof(volume_flux)) timestepper = $(typeof(timestepper)). ENSURE CORRECT BEFORE PROCEEDING")

# sleep(10)

result = ""

adaptive = true
N = 3
K = 64
abstol = 1e-6
reltol = 1e-4

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

include(MODULE)
result *= "$(sol.stats.naccept), $(sol.stats.nreject + sol.stats.naccept)\n"

N = 7
K = 32

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

include(MODULE)
result *= "$(sol.stats.naccept), $(sol.stats.nreject + sol.stats.naccept)\n"

abstol = 1e-8
reltol = 1e-6
N = 3
K = 64

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

include(MODULE)
result *= "$(sol.stats.naccept), $(sol.stats.nreject + sol.stats.naccept)\n"

N = 7
K = 32

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

include(MODULE)
result *= "$(sol.stats.naccept), $(sol.stats.nreject + sol.stats.naccept)\n"

abstol = 1e-10
reltol = 1e-8
N = 3
K = 64

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

include(MODULE)
result *= "$(sol.stats.naccept), $(sol.stats.nreject + sol.stats.naccept)\n"

N = 7
K = 32

if nameof(typeof(knapsack_solver)) == :ContinuousKnapsackSolver
    knapsack_solver = ContinuousKnapsackSolver((N + 2))
end

include(MODULE)
result *= "$(sol.stats.naccept), $(sol.stats.nreject + sol.stats.naccept)\n"

println(result)