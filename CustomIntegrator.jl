# Define a property for EEst with a custom setter
function Base.setproperty!(obj::OrdinaryDiffEq.ODEIntegrator, name::Symbol, value::Float64)
    if name === :EEst
        push!(total_error_estimates, value)
    end
    Core.setfield!(obj, name, value)
end

function Base.getproperty(obj::OrdinaryDiffEq.ODEIntegrator, name::Symbol)
    if preserve_positivity >= 0 && name === :dt
        global current_timestep = getfield(obj, name)
    end

    return getfield(obj, name)
end