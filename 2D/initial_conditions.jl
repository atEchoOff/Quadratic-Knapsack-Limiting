function initial_condition_density_wave(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    rho = 1 + .5 * sin(2 * pi * (x - .5 * t)) * sin(2 * pi * (y - t))
    v1 = .5
    v2 = 1.
    p = 1.0
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

function initial_condition_density_wave_low(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    rho = 1 + .98 * sin(2 * pi * (x - .5 * t)) * sin(2 * pi * (y - t))
    v1 = .5
    v2 = 1.
    p = 1.0
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

function initial_condition_khi(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    # KHI
    B = tanh(15 * y + 7.5) - tanh(15 * y - 7.5)
    rho = 0.5 + 0.75 * B
    u = 0.5 * (B - 1)
    v = 0.1 * sin(2 * pi * x)
    p = 1
    return prim2cons(SVector(rho, u, v, p), equations)
end

function initial_condition_football(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    if (0 <= x) & (0 <= y)
        rho = .5313
        u, v = 0, 0
        p = 0.4
    elseif (x < 0) & (0 <= y)
        rho = 1
        u, v = .7276, 0
        p = 1
    elseif (x < 0) & (y < 0)
        rho = .8
        u, v = 0, 0
        p = 1
    elseif (0 <= x) & (y < 0)        
        rho = 1
        u, v = 0, .7276
        p = 1
    else 
        @show x, y
    end

    return prim2cons(SVector(rho, u, v, p), equations)
end

function initial_condition_leaf(coords, t, ::CompressibleEulerEquations2D)
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

function initial_condition_riemann_clap(coords, t, ::CompressibleEulerEquations2D)
    x, y = coords

    if x >= 0 && y >= 0
        p = 1.1
        rho = 1.1
        u = 0.0
        v = 0.0
    elseif x < 0 && y >= 0
        p = .35
        rho = .01
        u = .8939
        v = 0.0
    elseif x < 0 && y < 0
        p = 1.1
        rho = 1.1
        u = .8939
        v = .8939
    else
        p = .35
        rho = .01
        u = 0.0
        v = .8939
    end

    return prim2cons(SVector(rho, u, v, p), equations)
end