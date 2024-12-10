include("L1_knapsack.jl")
include("L2_knapsack.jl")
include("L2_knapsack_weighted_a.jl")
L1 = ContinuousKnapsackSolver(2)
L2 = QuadraticKnapsackSolver{Float64}()
La = QuadraticKnapsackSolverA{Float64}()
b = .5

res = 1000
domain = range(-.2, 1, res)
values = Matrix{Float64}(undef, res, res)
for i in 1:res
    for j in 1:res
        x = domain[i]
        y = domain[j]

        aa = Float64[x, y]
        values[i, j] = La(aa, b)[2]
    end
end

heatmap(domain, domain, values)
xlabel!("a₁")
ylabel!("a₂")
