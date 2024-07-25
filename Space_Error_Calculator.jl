using LinearAlgebra
using Interpolations

function sort_and_unique(x::Vector{Float64}, y::Vector{Float64})
    # Combine x and y into a single matrix and sort by x values
    combined = sortperm(x)
    x_sorted = x[combined]
    y_sorted = y[combined]

    # Remove duplicates
    unique_indices = unique(i -> x[i], eachindex(x))
    x_unique = x_sorted[unique_indices]
    y_unique = y_sorted[unique_indices]

    return x_unique, y_unique
end

function estimate_L2_norm(x1::Matrix, f::Matrix, x2::Matrix, g::Matrix)::Float64
    # Flatten the matrices to vectors
    x1_flat = vec(x1)
    x2_flat = vec(x2)
    f_flat = vec(f)
    g_flat = vec(g)

    x1_flat, f_flat = sort_and_unique(x1_flat, f_flat)
    x2_flat, g_flat = sort_and_unique(x2_flat, g_flat)
    
    # Define the range of interpolation points
    common_points = range(-1 + 1/1000, stop=1 - 1/1000, length=1000)

    # Perform interpolation
    interp_f = LinearInterpolation(x1_flat, f_flat, extrapolation_bc=NaN)
    interp_g = LinearInterpolation(x2_flat, g_flat, extrapolation_bc=NaN)
    
    return 1 / 1000 * sum((interp_f.(common_points) - interp_g.(common_points)) .^2)
end