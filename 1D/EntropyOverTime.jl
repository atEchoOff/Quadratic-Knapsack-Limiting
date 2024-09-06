total_entropy(u) = sum(md.wJq .* entropy.(u, equations))

plot(sol.t, total_entropy.(sol.u))