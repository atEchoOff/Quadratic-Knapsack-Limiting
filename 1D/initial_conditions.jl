function initial_condition_linear_advection(x, t, equations::LinearScalarAdvectionEquation1D)
    return sin(pi * (x - t))
end

function initial_condition_burgers(x, t, equations::InviscidBurgersEquation1D)
    # u = exp(-10 * x^2)
    # if x < 0
    #     u = 0.0
    # else
    #     u = 1.0
    # end

    u = x

    return SVector(u)
end

function initial_condition_density_wave(x, t, equations::CompressibleEulerEquations1D)           
    # rho = 1.0 + exp(-100 * (x-0.25)^2)
    # # rho = 1.0 + .75 * (x > 0)
    # u = 0.0
    # p = rho^equations.gamma

    x = x - t
    x = x - floor(x)

    u = 1.
    rho = x + 3
    p = .1

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

function initial_condition_modified_sod(x, t, equations::CompressibleEulerEquations1D)
    if x[1] < .2
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