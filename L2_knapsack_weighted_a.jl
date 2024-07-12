struct QuadraticKnapsackSolverA{Ttol}
    tol::Ttol

    function QuadraticKnapsackSolverA{Ttol}(;tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol)
    end
end

function (s::QuadraticKnapsackSolverA)(a, b; upper_bounds=ones(length(a)))
    x = ones(length(a))
    return quadratic_knapsack_solverA!(x, a, b, upper_bounds; s.tol)
end

function quadratic_knapsack_solverA!(x, a, b, upper_bounds; tol = 100 * eps(), maxit=7) # maxit is a huge upper bound here. The iteration will take at most N - 1 iterations, but usually takes around 1 or 2
    I = upper_bounds .* (a .< tol)

    if a'I > b + tol
        # This is problem infeasibility
        # Here, we will just return the best we can do
        println("infeasibility (this should never happen in a well-posed problem)")
        return I
    end

    b = a' * upper_bounds - b

    if b < 0.0
        # x = 0 is the optimal solution for the minimization
        # meaning the optimal solution for the maximization problem is 1
        copy!(x, upper_bounds)
        return x
    end

    # Start the Newton iteration
    lambdak::Float64 = 0.0
    itercount = 0

    sgn_a = sign.(a)
    for _ in range(0, maxit)
        # Clip current solution within feasible domain
        x .= lambdak * sgn_a
        x = clamp.(x, 0.0, upper_bounds)

        f_val = a' * x - b

        if abs(f_val) < tol
            break
        end

        deriv = a' * (sgn_a .* (x .< upper_bounds) .* (lambdak * sgn_a .>= 0.0) .* (sgn_a .> 0.0))

        # Compute the next root
        lambdak -= f_val / deriv
        itercount += 1
    end

    x .= upper_bounds .- x

    # Comment these back in if you would like to test the iteration
    # @assert itercount <= 3
    # @assert a' * x <= original_b + tol
    # @assert abs(a' * x - original_b) <= tol
    # @assert all(x .<= upper_bounds)
    # @assert all(x .>= 0.0)
    return x
end