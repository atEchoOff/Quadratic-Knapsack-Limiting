MODULE = "burgers.jl"

include("Running_Interface.jl")
Li = Float64[]
Ki = Float64[]

for i in 2 .^ ((40:70) * .1)
    global K
    K = floor(Int, i)
    push!(Ki, 1 / K)

    include(MODULE)

    nm = L2_error
    println("Obtained norm $nm")
    push!(Li, nm)
end

println(fit(log.(Ki), log.(Li), 1))