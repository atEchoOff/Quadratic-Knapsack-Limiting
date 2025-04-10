using LinearAlgebra
using ForwardDiff

function roll_up(u, m, n, d)
    return reshape(reinterpret(SVector{d, Float64}, u), m, n)
end

function unroll(u)
    return vec(reinterpret(Float64, u))
end

function get_jacobian_rhs(du, u; epsilon=1e-6)
    m = size(u, 1)
    n = size(u, 2)
    d = length(u[1, 1])

    A = Matrix{Float64}(undef, m * n * d, m * n * d)

    for i in 1:size(A, 2)
        ei = 1. * (1:(m * n * d) .== i)
        u_forward = u + epsilon * roll_up(ei, m, n, d)
        u_backward = u - epsilon * roll_up(ei, m, n, d)

        du_forward = copy(unroll(rhs!(du, u_forward, cache, 0.)))
        du_backward = copy(unroll(rhs!(du, u_backward, cache, 0.)))

        A[:, i] .= 1 / (2 * epsilon) * unroll(du_forward - du_backward)
    end

    return A
end

function get_jesse_jacobian_rhs(du, u)
    rhs_fd = let cache=cache
        function rhs_fd(u_in::AbstractArray{T}) where {T}
            u = reinterpret(reshape, SVector{nvariables(equations), T}, 
                            reshape(u_in, :, size(md.x)...))
            du = similar(u)
            rhs!(du, u, cache, 0.0)
            return reinterpret(reshape, T, du)
        end
    end

    u_in = reinterpret(reshape, Float64, u)
    return ForwardDiff.jacobian(rhs_fd, u_in)
end