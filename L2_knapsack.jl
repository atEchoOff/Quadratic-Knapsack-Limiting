struct QuadraticKnapsackSolver{Ttol}
    tol::Ttol

    function QuadraticKnapsackSolver{Ttol}(;tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol)
    end
end

function (s::QuadraticKnapsackSolver)(a, b; upper_bounds=ones(length(a)), w=ones(length(a)))
    # Take a, b optimization parameters
    # Take upper_bounds (default is ones)
    # Take w weights (default is ones)
    # Return (not by parameter) optimal θ
    x = ones(length(a))
    return quadratic_knapsack_solver!(x, a, b, upper_bounds, w; s.tol)
end

function smooth_clamp(x, b)
    x = x / b
    if x < 0
        return 0
    elseif x > 1
        return 1
    elseif x < .1
        return -100x^3 + 20x^2
    elseif x > .9
        return 1 - (-100*(1-x)^3 + 20*(1-x)^2)
    else
        return x
    end
end

function quadratic_knapsack_solver!(x, a, b, upper_bounds, w; tol = 100 * eps(), maxit=200)
    # maxit is a huge upper bound here. The iteration will take at most N + 1 iterations, but usually takes around 1 - 3, and rarely 4.
    # Note also, if maxit is reached, it likely implies that b is negative, so the problem needs cleaning

    ### Optional infeasibility check, uncomment to see if issues are caused by infeasibility
    # I = upper_bounds .* (a .< tol)

    # if (a'I > b + tol) || (b < 0)
    #     # This is problem infeasibility
    #     # Here, we will just return the best we can do
    #     println("infeasibility (this should never happen in a well-posed problem)")
    #     return I
    # end

    # Problem infeasibility... this is bad
    # Just return low order solution... best we can do
    if b <= 0.0
        x .= 0.0
        # println("Infeasibility Detected")
        return x
    end

    # original_b = b # for sanity check
    b = a' * upper_bounds - b

    if b <= 0.0
        # x = 0 is the optimal solution for the minimization
        # meaning the optimal solution for the maximization problem is upper_bounds
        copy!(x, upper_bounds)
        return x
    end

    # Start the Newton iteration
    lambdak::Float64 = 0.0
    itercount = 0 # for sanity check

    a_over_w = a ./ w

    # sgn_a = sign.(a_over_w)
    # a_over_w .= abs.(a_over_w)
    # a .= abs.(a)
    
    for _ in range(0, maxit)
        # Clip current solution within feasible domain
        x .= lambdak * a_over_w
        x = clamp.(x, 0.0, upper_bounds)

        f_val = a' * x - b

        if abs(f_val) < tol
            break
        end

        deriv = a' * (a_over_w .* (x .< upper_bounds) .* (lambdak * a_over_w .>= 0.0) .* (a_over_w .> 0.0))

        # Compute the next root
        lambdak -= f_val / deriv

        itercount += 1 # for sanity check
    end

    # x .= sgn_a .* x

    x .= upper_bounds .- x

    ### These are all my sanity checks. Non well-posed problems may break them, so if issues are found, uncomment these and the sanity check comments above for checking.
    if itercount > 3
        println("The itercount was $itercount")
    end

    # @assert all(x .== 1.0)

    # @assert itercount <= 4

    # @assert a' * x <= original_b + tol
    # @assert abs(a' * x - original_b) <= tol
    # @assert all(x .<= upper_bounds)
    # @assert all(x .>= 0.0)
    return x
end