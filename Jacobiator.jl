using LinearAlgebra

function roll_up(u, m, n, d)
    return reshape(reinterpret(SVector{d, Float64}, u), m, n)
end

function unroll(u)
    return vec(reinterpret(Float64, u))
end

function get_jacobian_rhs(du, u; epsilon=1e-3)
    m = size(u, 1)
    n = size(u, 2)
    d = length(u[1, 1])

    A = Matrix{Float64}(undef, m * n * d, m * n * d)

    for i in 1:size(A, 2)
        ei = 1. * (1:(m * n * d) .== i)
        u_forward = u + epsilon * roll_up(ei, m, n, d)
        u_backward = u - epsilon * roll_up(ei, m, n, d)



        A[:, i] .= 1 / (2 * epsilon) * unroll(rhs!(du, u_forward, cache, 0.) - rhs!(du, u_backward, cache, 0.))
    end

    return A
end