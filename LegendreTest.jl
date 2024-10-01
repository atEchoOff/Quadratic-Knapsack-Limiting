using LegendrePolynomials

function lpol(deg::Integer)
    return function inner(x)
        return Pl(x, deg)
    end
end

V = Matrix{Float64}(undef, 8, 8)

for i in 1:8
    if i == 1
        V[:, i] .= 0
    else
        V[:, i] = lpol(i - 1).(rd.rq)
    end
end