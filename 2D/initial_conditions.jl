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

function initial_condition_khi_asymmetric(coords, t, equations::CompressibleEulerEquations2D)
    x, y = coords

    # A = .8
    A = 3 / 7
    rho1 = 0.5 * one(A) # recover original with A = 3/7
    rho2 = rho1 * (1 + A) / (1 - A)

    # B is a discontinuous function with value 1 for -.5 <= x <= .5 and 0 elsewhere
    slope = 15
    B = 0.5 * (tanh(slope * y + 7.5) - tanh(slope * y - 7.5))

    rho = rho1 + B * (rho2 - rho1)  # rho ∈ [rho_1, rho_2]
    v1 = B - 0.5                    # v1  ∈ [-.5, .5]
    # v2 = 0.1 * sin(2 * pi * x[1]) 
    v2 = 0.1 * sin(2 * pi * x) * (1 + .01 * sin(pi * x) * sin(pi * y)) # symmetry breaking
    p = 1.0
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

function initial_condition_khi_dissipated(coords, t, equations::CompressibleEulerEquations2D)
    x, y = coords

    # A = .8
    A = 3 / 7
    rho1 = 0.5 * one(A) # recover original with A = 3/7
    rho2 = rho1 * (1 + A) / (1 - A)

    # B is a discontinuous function with value 1 for -.5 <= x <= .5 and 0 elsewhere
    slope = 15
    B = 0.5 * (tanh(slope * y + 7.5) - tanh(slope * y - 7.5))

    rho = rho1 + B * (rho2 - rho1)  # rho ∈ [rho_1, rho_2]
    v1 = B - 0.5                    # v1  ∈ [-.5, .5]
    # v2 = 0.1 * sin(2 * pi * x[1]) 
    v2 = 0.1 * sin(2 * pi * x)
    p = 1.0
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

function initial_condition_sedov_blastwave(coords, t, equations::CompressibleEulerEquations2D)
    x, y = coords
    h = 3 / K # mesh size
    r0 = 4h

    r = sqrt(x^2 + y^2)
    if r < r0
        rho = 1
        v1 = 0
        v2 = 0
        p = .4 / (pi * r0^2)
    else
        rho = 1
        v1 = 0
        v2 = 0
        p = 1e-5
    end

    return prim2cons(SVector(rho, v1, v2, p), equations)
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