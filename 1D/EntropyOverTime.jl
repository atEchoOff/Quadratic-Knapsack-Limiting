total_entropy(u) = sum(md.wJq .* Trixi.entropy.(u, equations))

plot(sol.t, total_entropy.(sol.u))