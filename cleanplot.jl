function dectoidx(dec, K)
    return floor(Int, round((dec + 1) * 5K / 2) + 1)
end

function plot2D(x, y, u)
    @assert length(x) == length(y)
    @assert length(y) == length(u)

    x_idx = dectoidx.(x, K)
    y_idx = dectoidx.(y, K)
    
    x_max = maximum(x_idx)
    y_max = maximum(y_idx)

    @assert x_max == y_max

    u_plot = zeros(Float64, x_max, y_max)

    for i in range(1, length(u))
        u_plot[x_idx[i], y_idx[i]] = u[i]
    end

    heatmap(u_plot)
end