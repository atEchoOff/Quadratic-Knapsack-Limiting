struct QuadraticKnapsackSolverA{Ttol}
    tol::Ttol

    function QuadraticKnapsackSolverA{Ttol}(;tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol)
    end
end

function (s::QuadraticKnapsackSolverA)(a, b; upper_bounds=ones(length(a)))
    x = ones(length(a))
    return quadratic_knapsack_solverA!(x, a, b, upper_bounds)
end

function quadratic_knapsack_solverA!(x, a, b, upper_bounds)
    # Problem infeasibility... this is bad
    # Just return low order solution... best we can do
    if b <= 0.0
        x .= 0.0
        # println("Infeasibility Detected")
        return x
    end

    if a'upper_bounds <= b
        x .= upper_bounds
        return x
    end

    atupper_bounds = a'upper_bounds
    # sgn_a = sign.(a)
    # a .= abs.(a)

    lambda = (atupper_bounds - b) / (upper_bounds' * (a .* (a .> 0.0)))

    x .= upper_bounds .* clamp.(lambda * sign.(a), 0.0, 1.0)

    # x .= sgn_a .* x

    x .= upper_bounds - x
    return x
end