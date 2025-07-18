CFLs = []
Δx = .02
for i in 1:(length(sol.t) - 1)
    ui = cons2prim.(sol.u[i], equations)
    rho = getindex.(ui, 1)
    v1 = getindex.(ui, 2)
    p = getindex.(ui, 3)
    speeds = @. abs(v1) + sqrt(1.4 * p / rho)
    push!(CFLs, 4 * dt / Δx * maximum(speeds ./ rd.wq))
end