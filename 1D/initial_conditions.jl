using Distributions
using SodShockTube

# Stop being stupid Julia
solve = OrdinaryDiffEq.solve

function initial_condition_linear_advection(x, t, equations::LinearScalarAdvectionEquation1D)
    # return sin(pi * (x - t))
    return x
end

function initial_condition_burgers(x, t, equations::InviscidBurgersEquation1D)
    u = exp(-10 * x^2)
    # u = sin(π * x - .7) + 2
    # if x < 0
    #     u = .5
    # else
    #     u = 1.5
    # end

    # u = sign(x)

    return SVector(u)
end

function initial_condition_density_wave(x, t, equations::CompressibleEulerEquations1D; dist=Normal(0., 0.))           
    # rho = 1.0 + exp(-100 * (x-0.25)^2)
    # # rho = 1.0 + .75 * (x > 0)
    # u = 0.0
    # p = rho^equations.gamma

    # x = x - t
    # x = x - floor(x)

    u = 1.
    # rho = sin(pi * (x - t))^2 + .1
    rho = sin(8pi * (x - t))^2 + .1
    # rho = 1 + .98 * sin(2π * x)
    # rho = sin(pi * (x - t))^2 + .1 + rand(dist)
    # rho = (x - t) - floor(x - t) + .1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_modified_sod(x, t, equations::CompressibleEulerEquations1D)

    if t > 0.
        # Instead, return the analytical solution (thanks to SodShockTube repo)
        # Set up a shock tube problem
        # This is set up for Sod Shock, not modified
        problem = ShockTubeProblem(
            geometry = (0.0, 1.0, .5), # left edge, right edge, initial shock location
            left_state = (ρ = 1.0, u = 0., p = 1.0),
            right_state = (ρ = 0.125, u = 0.0, p = 0.1),
            t = t,
            γ = 1.4
        );

        _, _, values = SodShockTube.solve(problem, [x])
        return prim2cons(SVector{3, Float64}(values.ρ[1], values.u[1], values.p[1]), equations)
    end
    if x[1] < .3
        rho = 1.0
        v1 = .75
        p = 1.0
    else
        rho = .125
        v1 = 0.0
        p = .1
    end
    return prim2cons(SVector(rho, v1, p), equations)
end