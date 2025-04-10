using LsqFit
using LaTeXStrings

function fit_to_order(x_data, y_data, order)
    @. model(x, p) = p[1] * x ^ order

    return curve_fit(model, x_data, y_data, [.5]).param[1]
end

scatter(h1, E1, xaxis=:log, yaxis=:log, label = L"$N$ = 1")
scatter!(h2, E2, label = L"$N$ = 2")
scatter!(h3, E3, label = L"$N$ = 3")
scatter!(h4, E4, label = L"$N$ = 4")

plot!(legend=:topleft)

xlabel!(L"h")
ylabel!(L"$L^2$ Error")

o1 = fit_to_order(h1, E1, 2)
o2 = fit_to_order(h2, E2, 3)
o3 = fit_to_order(h3, E3, 4)
o4 = fit_to_order(h4, E4, 5)

plot!(annotations=([2e-2], [1e-3], text(L"$\mathcal{O}\left(h^2\right)$", :red)))
plot!(h1, o1 * (h1 .^ 2), label = "", c=:red, lw = 2)

plot!(annotations=([5e-2], [2e-5], text(L"$\mathcal{O}\left(h^3\right)$", :red)))
plot!(h2, o2 * (h2 .^ 3), label = "", c=:red, lw = 2)

plot!(annotations=([1.4e-1], [5e-6], text(L"$\mathcal{O}\left(h^4\right)$", :red)))
plot!(h3, o3 * (h3 .^ 4), label = "", c=:red, lw = 2)

plot!(annotations=([2.3e-1], [1e-6], text(L"$\mathcal{O}\left(h^5\right)$", :red)))
plot!(h4, o4 * (h4 .^ 5), label = "", c=:red, lw = 2)

niceplot!()