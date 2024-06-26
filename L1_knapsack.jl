using LinearAlgebra

struct ContinuousKnapsackSolver{T, Ttol}
    indices_sorted::T # tmp storage for sorting
    tol::Ttol
end

ContinuousKnapsackSolver(a::AbstractVector; tol = 100 * eps()) = 
    ContinuousKnapsackSolver(length(a); tol)

function ContinuousKnapsackSolver(n::Int; tol=100 * eps()) 
    indices_sorted = [collect(1:n) for _ in 1:Threads.nthreads()]
    return ContinuousKnapsackSolver(indices_sorted, tol)
end

function (s::ContinuousKnapsackSolver)(x, a, b)
    return knapsack_solver!(x, s.indices_sorted[Threads.threadid()], a, b; s.tol)
end

function (s::ContinuousKnapsackSolver)(a, b; upper_bounds=ones(length(a)))
    x = copy(upper_bounds)
    return knapsack_solver!(x, s.indices_sorted[Threads.threadid()], a, b; s.tol)
end

# solves the problem
# max ∑ x_i
# s.t. a' * x <= b : ∑ⁿ a_j * x_j = O(n)
#      0 ≤ x_i ≤ 1 (or in general, (x_max)_i)
# using an exact greedy algorithm 
function knapsack_solver!(x, indices_sorted, a, b; tol = 100 * eps())
    if length(a) < 20
        sortperm!(indices_sorted, a, alg=InsertionSort, rev=true) 
    else
        sortperm!(indices_sorted, a, alg=QuickSort, rev=true)        
    end

    upper_bounds = copy(x)

    s = dot(a, x)
    for i in indices_sorted
        # ignore negative a[i] since we can take x[i] = x_i_max = 1.0 
        if a[i] > tol
            # remove next largest positive contribution
            s = s - a[i] * x[i] 
            if s > b + tol || a[i] < tol
                # if s without x[i] is still > b, remove the contribution by setting x[i] = 0            
                x[i] = 0.0
            else
                # if s without x[i] ≤ b, add the max that we can back in
                x[i] = min(upper_bounds[i], (b - s) / a[i]) # want to find x[i] such that s + a[i]* x[i] = b
            end
        end        
    end

    return x
end


# using JuMP
# import HiGHS


# # use JuMP to compute a reference solution
# function solve_knapsack_problem(;
#     weight::Vector{T},
#     capacity::T,
# ) where {T<:Real}
#     N = length(weight)

#     model = Model(HiGHS.Optimizer)
#     set_silent(model)

#     # Declare the decision variables as between 0, 1
#     @variable(model, 0 <= x[1:N] <= 1)

#     # Objective: maximize profit.
#     @objective(model, Max, sum(x))

#     # Constraint: can carry all items.
#     @constraint(model, weight' * x <= capacity)

#     # Solve problem using MIP solver
#     optimize!(model)
#     return value.(x)
# end


# # test knapsack solvers
# a = rand(50) + 0.5 * randn(50)
# b = rand()

# x_jump = solve_knapsack_problem(; weight = a, capacity = b)

# knapsack_solve! = ContinuousKnapsackSolver(a)
# knapsack_solve!(x, a, b)
# # x = knapsack_solver(a, b)

# @show norm(x - x_jump)