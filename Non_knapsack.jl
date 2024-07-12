struct NonKnapsackSolver{Ttol}
    tol::Ttol

    function NonKnapsackSolver{Ttol}(;tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol)
    end
end

function (s::NonKnapsackSolver)(a, b; upper_bounds=ones(length(a)), w=ones(length(a)))
    # Take a, b optimization parameters
    # Take upper_bounds (default is ones)
    # Take w weights (default is ones)
    # Return (not by parameter) optimal θ
    x = ones(length(a))
    return non_knapsack_solver!(x, a, b, upper_bounds, w; s.tol)
end

function non_knapsack_solver!(x, a, b, upper_bounds, w; tol = 100 * eps(), maxit=200)
    return upper_bounds
end