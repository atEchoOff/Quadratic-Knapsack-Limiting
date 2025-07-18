include("L1_knapsack.jl")
include("L2_knapsack_minimizer.jl")
include("L2_knapsack_weighted_a.jl")
L1 = ContinuousKnapsackSolver(2)
L2 = QuadraticKnapsackMinimizer{Float64}()
La = QuadraticKnapsackSolverA{Float64}()
b = .5

epsilon = 1e-6
e1 = epsilon * [1, 0]
e2 = epsilon * [0, 1]

res = 1000
domain = range(.25, 1.2, res)
values = Matrix{Float64}(undef, res, res)
for i in 1:res
    for j in 1:res
        x = domain[i]
        y = domain[j]

        aa = Float64[y, x]
        values[i, j] = L2(aa, b)[2]
        # values[i, j] = ((L2(aa + e2, b) - L2(aa - e2, b)) / (2epsilon))[2]
    end
end

heatmap(domain, domain, values)
xlabel!(L"$a_1$")
ylabel!(L"$a_2$")
niceplot!()