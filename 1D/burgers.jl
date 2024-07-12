using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots

include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_weighted_a.jl")

function rhs!(du, u, cache, t)
    (; physics, rd, md) = cache
    (; Q_skew) = cache.high_order_operators
    (; Q_skew_low) = cache.low_order_operators    
    (; Δ, R) = cache.fv_operators
    (; equations, volume_flux, surface_flux) = physics
    (; VDM_inv, shock_capturing, nodewise_shock_capturing) = cache

    # interface flux
    uf = rd.Vf * u
    uP = uf[md.mapP]

    interface_flux = @. surface_flux(uf, uP, SVector{1}(md.nxJ), equations)
    du .= rd.LIFT * interface_flux
   
    rhs_vol_high = similar(du, size(du, 1))
    rhs_vol_low = similar(rhs_vol_high)
    
    A = fill!(similar(du, size(du, 1), size(du, 1)), zero(eltype(du)))
    a = zeros(size(du, 1), size(du, 1))
    for e in 1:md.num_elements

        fill!(rhs_vol_high, zero(eltype(rhs_vol_high)))
        fill!(rhs_vol_low, zero(eltype(rhs_vol_low)))
    
        # high order update
        u_element = view(u, :, e)
        v = cons2entropy.(u_element, equations)
        for i in eachindex(u_element), j in eachindex(u_element)            
            if i > j
                fij = volume_flux(u_element[i], u_element[j], 1, equations)
                FH_ij = Q_skew[i, j] * fij
                rhs_vol_high[i] += FH_ij
                rhs_vol_high[j] += -FH_ij

                FL_ij = zero(eltype(rhs_vol_low))
                if abs(Q_skew_low[i, j]) > 1e-12
                    nij = Q_skew_low[i, j] / abs(Q_skew_low[i, j])
                    fij = surface_flux(u_element[i], u_element[j], SVector{1}(nij), equations)
                    FL_ij = abs(Q_skew_low[i, j]) * fij
                    rhs_vol_low[i] += FL_ij
                    rhs_vol_low[j] -= FL_ij    
                end

                # create antidiffusive flux matrix
                A[i,j] = FH_ij - FL_ij
                A[j,i] = -A[i, j]
            end
        end

        (; blend) = cache
        
        # optimization targets
        # println(cache.B)
        # println(psi.(u_element, equations))
        # println(cache.B * psi.(u_element, equations))
        # println(sum(dot.(v, rhs_vol_low)))
        b = sum(cache.B * psi.(u_element, equations))[1] + sum(dot.(v, rhs_vol_low))

        ### Shock capturing optimization target
        if shock_capturing
            μ = VDM_inv * getindex.(u_element, 1)
            sum_to_end = sum(μ.^2)

            s = log10(max(μ[end]^2 / sum_to_end, μ[end - 1]^2 / (sum_to_end - μ[end]^2)))
            κ = 1.0
            s0 = -4 * log10(size(u_element, 1) - 1)
            ϵ = 0.0
            if s < s0 + κ
                ϵ = 1.0
            elseif s < s0 - κ
                ϵ = 0.0
            else
                ϵ = 0.5 * (1.0 - sin(π * (s - s0) / (2 * κ)))
            end

            b = (1 - .01 * ϵ) * b
        end

        if blend == :elementwise
            # 1. elementwise blending of rhs_vol_high/low
            denom = (sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(1.0, θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)

        elseif blend == :subcell
            # 2. subcell Δ, R blending of rhs_vol_high/low        
            # @assert b < 100 * eps()
            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))
            a = vec(a)

            θ = cache.knapsack_solver(a, b)

            if nodewise_shock_capturing
                epsilon = 120.0
                # @. a *= -log(1/e * θ)
                # @. a *= 1 / cbrt(θ)
                @. a *= 1 - 4 * N^2 * log(θ)
                # @. a *= 1 + N*sqrt(-log(θ))
                # # @. a *= (1 - 2 * N * (K / 64) * log(θ + 1e-8))
                # @. a *= -epsilon * θ + (1 + epsilon)
                # # @. a *= (-sqrt(epsilon) * θ + sqrt(epsilon))^2 + 1
                # @. a *= -epsilon * min(1, θ + 1e-1) + (1 + epsilon)

                θ = cache.knapsack_solver(a, b)
                # @. a *= -epsilon * θ + (1 + epsilon)
                # # @. a *= (-sqrt(epsilon) * θ + sqrt(epsilon))^2 + 1

                # θ = cache.knapsack_solver(a, b)
            end

            view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ) * R * rhs_vol_high + Diagonal(1 .- θ) * R * rhs_vol_low))

        else

            view(du, :, e) .+= rd.M \ rhs_vol_high
        end
    end
    
    # invert Jacobian and mass matrix
    du ./= -md.J
end

psi(u, ::InviscidBurgersEquation1D) = 1/6 * u .^ 3

function initial_condition_basic(x, t, equations::InviscidBurgersEquation1D)
    u = exp(-10 * x^2)

    return SVector(u)
end

initial_condition = initial_condition_basic

initial_limiting_coefficients(R, u; delta=0) = 1 .- delta * rand(size(R, 1))

rd = RefElemData(Line(), SBP(), N) 
(VX, ), EToV = uniform_mesh(Line(), K)

md = MeshData((VX,), EToV, rd)
md = make_periodic(md)

equations = InviscidBurgersEquation1D()

u0 = rd.Pq * initial_condition.(md.xq, 0.0, equations)
u = copy(u0)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

# subcell FV operators
Δ = diagm(0 => -ones(rd.N+1), 1=>ones(rd.N+1))[1:end-1,:]
R = tril(ones(size(Δ')), -1)

cache = (; 
           rd, md, 
           B = Diagonal([-1; zeros(rd.N-1); 1]),
           high_order_operators=(; Q_skew = rd.M * rd.Dr - (rd.M * rd.Dr)'), 
           low_order_operators=(; Q_skew_low = Qr-Qr'), 
           fv_operators = (; Δ, R), 
           blend = blend,
           knapsack_solver = knapsack_solver,
        physics = (; equations, volume_flux = volume_flux, surface_flux = flux_lax_friedrichs),
        shock_capturing = shock_capturing,
        nodewise_shock_capturing = nodewise_shock_capturing,
        VDM_inv = inv(rd.VDM)
        )

# dt = 0.25 * estimate_h(rd, md) * (1 / (2 * rd.N + 1))
tspan = (0, 3.0)
ode = ODEProblem(rhs!, u, tspan, cache)
@time sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, saveat=LinRange(tspan..., 10), callback=AliveCallback(alive_interval=1000), adaptive=adaptive)

# @gif for u in sol.u
#     plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0, 2))
# end
u = sol.u[end]
# plot!(rd.Vp * md.x, getindex.(initial_condition.(rd.Vp * md.x, tspan[2], equations), 1), leg=false)

rq, wq = gauss_quad(0, 0, N+2)
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x

L1_error = sum(wJq .* map(x -> sum(abs.(x)), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end]))
L2_error = sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end])))

println("N = $N, K1D = $(md.num_elements), L1_error = $L1_error, L2_error = $L2_error")


plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false)