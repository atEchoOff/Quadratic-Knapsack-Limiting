function min_norm_solution(f, b, n)
    # Create the matrix representation of f by applying f to the standard basis vectors
    A = zeros(eltype(b), length(b), n)
    for i in 1:n
        e = zeros(n)
        e[i] = 1.0
        A[:, i] = f(e)
    end

    return A \ b
end