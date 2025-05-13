CFLs = []
Δx = .03 # For Sedov Blast
for i in 1:(length(sol.t) - 1)
    ui = cons2prim.(sol.u[i], equations)
    rho = getindex.(ui, 1)
    v1 = getindex.(ui, 2)
    v2 = getindex.(ui, 3)
    p = getindex.(ui, 4)
    speeds = @. norm(zip(v1, v2)) + sqrt(1.4 * p / rho)
    push!(CFLs, 4 * dt / Δx * maximum(speeds ./ rd.wq))
end