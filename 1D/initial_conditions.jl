using SodShockTube

# Stop being stupid Julia
solve = OrdinaryDiffEq.solve

function initial_condition_linear_advection(x, t, equations::LinearScalarAdvectionEquation1D)
    return exp(-10 * x^2)
end

function initial_condition_burgers(x, t, equations::InviscidBurgersEquation1D)
    u = exp(-10 * x^2)
    return SVector(u)
end

function initial_condition_density_wave_jesse(x, t, equations::CompressibleEulerEquations1D)
    u = .1
    rho = .98 * sin(2pi * (x - .1t)) + 1
    p = 10.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_density_wave_slow(x, t, equations::CompressibleEulerEquations1D)
    u = .1
    rho = .98 * sin(pi * (x - .1t)) + 1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_density_wave(x, t, equations::CompressibleEulerEquations1D)
    u = 1.
    rho = .5 * sin(2 * pi * (x - t)) + 1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_density_wave_fast(x, t, equations::CompressibleEulerEquations1D)
    u = 1.7
    rho = .5 * sin(pi * (x - 1.7t)) + 1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_density_wave_low(x, t, equations::CompressibleEulerEquations1D)
    u = 1.
    rho = .98 * sin(2 * pi * (x - t)) + 1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_constant_one(x, t, equations::CompressibleEulerEquations1D)
    u = .1
    rho = 1.
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_sod_shock(x, t, equations::CompressibleEulerEquations1D)
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
    else
        if x[1] < .5
            rho = 1.0
            v1 = 0.
            p = 1.0
        else
            rho = .125
            v1 = 0.0
            p = .1
        end
    end

    return prim2cons(SVector(rho, v1, p), equations)
end

function initial_condition_leblanc_shocktube(x, t, equations::CompressibleEulerEquations1D)
    if x[1] <= 0
        rho = 2.
        v1 = 0.
        p = 1e9
    else
        rho = .001
        v1 = 0.
        p = 1.
    end
    return prim2cons(SVector(rho, v1, p), equations)
end

function initial_condition_modified_sod(x, t, equations::CompressibleEulerEquations1D)
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

function initial_condition_modified_squared_sod(x, t, equations::CompressibleEulerEquations1D)
    if x[1] < .3
        rho = 1.0
        v1 = .75
        p = 1.0
    else
        rho = .0125
        v1 = 0.0
        p = .01
    end
    return prim2cons(SVector(rho, v1, p), equations)
end

function initial_condition_shu_osher(x, t, equations::CompressibleEulerEquations1D)
    if x[1] < -4
        rho = 3.857143
        v1 = 2.629369
        p = 10.3333 
    else
        rho = 1 + .2 * sin(5 * x[1])
        v1 = 0.0
        p = 1.0
    end
    return prim2cons(SVector(rho, v1, p), equations)
end