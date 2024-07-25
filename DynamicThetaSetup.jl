# Define a property for EEst with a custom setter
# function Base.setproperty!(obj::OrdinaryDiffEq.ODEIntegrator, name::Symbol, value::Float64)
#     if name === :EEst
#         push!(total_error_estimates, value)
#     end
#     Core.setfield!(obj, name, value)
# end

function Base.:*(a::SVector, b::SVector)
    return a .* b
end

function vector_to_svector_vector(flat_vector::Vector{T}, dim::Int) where T
    n = length(flat_vector)
    if n % dim != 0
        throw(ArgumentError("The length of the flat vector must be a multiple of the specified dimension"))
    end
    num_sv = n ÷ dim
    svectors = Vector{SVector{dim, T}}(undef, num_sv)
    for i in 1:num_sv
        svectors[i] = SVector{dim, T}(flat_vector[(i-1)*dim+1:i*dim]...)
    end
    return svectors
end