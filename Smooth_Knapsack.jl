using Roots

struct SmoothKnapsackSolver{Ttol}
    tol::Ttol

    function SmoothKnapsackSolver{Ttol}(;tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol)
    end
end

function (s::SmoothKnapsackSolver)(a, b; upper_bounds=ones(length(a)), w=ones(length(a)))
    # Take a, b optimization parameters
    # Take upper_bounds (default is ones)
    # Take w weights (default is ones)
    # Return (not by parameter) optimal θ
    x = ones(length(a))
    return smooth_knapsack_solver!(x, a, b, upper_bounds, w; s.tol)
end

function sigmoid(x::Real)
    return 1 / (1 + exp(-x))
end

function smooth_knapsack_solver!(x, a, b, upper_bounds, w; tol = 100 * eps(), maxit=200)
    if sum(a) <= b + tol
        x .= 1.0
        return x
    end
    
    f(c) = a' * sigmoid.(-2 * (a .- c)) - b

    c = find_zero(f, (-10, 10), Roots.Brent())

    x .= sigmoid.(-2 * (a .- c))
    return x
end