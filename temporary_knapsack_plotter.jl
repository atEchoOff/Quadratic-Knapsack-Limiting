include("L1_knapsack.jl")
include("L2_knapsack_minimizer.jl")
include("L2_knapsack_weighted_a.jl")
L1 = ContinuousKnapsackSolver(2)
L2 = QuadraticKnapsackMinimizer{Float64}()
La = QuadraticKnapsackSolverA{Float64}()

epsilon = 1e-6
e1 = epsilon * [1, 0]
e2 = epsilon * [0, 1]

res = 1000
domain = range(.25, 1.2, res)
brange = range(-.1, .5, 100)

values = Array{Float64}(undef, res, res, res)
for i in 1:res
    for j in 1:res
        for k in 1:100
            x = domain[i]
            y = domain[j]
            b = brange[k]

            aa = Float64[x, y]
            values[i, j, k] = L1(aa, b)[2]
            # values[i, j] = ((L2(aa + e2, b) - L2(aa - e2, b)) / (2epsilon))[2]
        end
    end
end

@gif for k in 1:100
    heatmap(domain, domain, values[:, :, k], clims=(0, 1))
    xlabel!(L"$a_1$")
    ylabel!(L"$a_2$")
    title!(L"$b$ = " * string(round(brange[k], sigdigits=2)))
    niceplot!()
end