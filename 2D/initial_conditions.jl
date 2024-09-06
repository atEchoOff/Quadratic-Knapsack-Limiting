function initial_condition_density_wave(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    rho = 1 + sin(pi * (x - .5 * t))^2 * sin(pi * (y - t))^2
    v1 = .5
    v2 = 1.
    p = 1.0
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

function initial_condition_khi(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    # rho = 1 + 0.5 * (abs(x) < 0.5) * (abs(y) < 0.5)
    # v1 = 0.
    # v2 = 0.
    # p = rho^equations.gamma
    # return prim2cons(SVector(rho, v1, v2, p), equations)

    # KHI
    B = tanh(15 * y + 7.5) - tanh(15 * y - 7.5)
    rho = 0.5 + 0.75 * B
    u = 0.5 * (B - 1)
    v = 0.1 * sin(2 * pi * x)
    p = 1
    return prim2cons(SVector(rho, u, v, p), equations)
end

function initial_condition_riemann(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    if x >= 0 && y >= 0
        p = 1.1
        rho = 1.1
        u = 0.0
        v = 0.0
    elseif x < 0 && y >= 0
        p = .35
        rho = .5065
        u = .8939
        v = 0.0
    elseif x < 0 && y < 0
        p = 1.1
        rho = 1.1
        u = .8939
        v = .8939
    else
        p = .35
        rho = .5065
        u = 0.0
        v = .8939
    end

    return prim2cons(SVector(rho, u, v, p), equations)
end