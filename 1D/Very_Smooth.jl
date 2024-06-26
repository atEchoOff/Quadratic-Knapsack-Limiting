using StartUpDG
using StaticArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots

include("../L2_knapsack.jl")

function rhs!(du, u, cache, t)
    (; physics, rd, md) = cache
    (; Q_skew) = cache.high_order_operators
    (; Q_skew_low) = cache.low_order_operators    
    (; Δ, R) = cache.fv_operators
    (; equations, volume_flux, surface_flux) = physics

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
        
        # optimization target
        b = sum(cache.B * psi.(u_element, equations)) + sum(dot.(v, rhs_vol_low))

        # initial limiting coefficients
        θ_init = initial_limiting_coefficients(R, u)

        if blend == :elementwise
            # 1. elementwise blending of rhs_vol_high/low
            denom = (sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(minimum(θ_init), θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)

        elseif blend == :subcell
            # 2. subcell Δ, R blending of rhs_vol_high/low
            a = v' * Δ * Diagonal(R * (rhs_vol_low - rhs_vol_high))

            # Call knapsack
            θ = cache.knapsack_solver(vec(a), b; upper_bounds=θ_init)
            view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ) * R * rhs_vol_high + 
                                            Diagonal(1 .- θ) * R * rhs_vol_low))
        else
            view(du, :, e) .+= rd.M \ rhs_vol_high
        end
    end
    
    # invert Jacobian and mass matrix
    du ./= -md.J
end

psi(u, ::CompressibleEulerEquations1D) = u[2] # rho * v1

function initial_condition_basic(x, t, equations::CompressibleEulerEquations1D)           
    # rho = 1.0 + exp(-100 * (x-0.25)^2)
    # rho = 1.0 + 1.0 * (x > 0)
    # u = 0.0
    # p = rho^equations.gamma 

    u = 1.0
    rho = 1 + .5 * sin(2 * pi * (x - u * t))
    p = .1

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

# Number of nodes in each element
N = 8
blend = :subcell

# Set some arbitrary limiting coefficients. These are just to test stability of the solver when we have strict upper bounds
initial_limiting_coefficients(R, u; delta=1e-2) = 1 .- delta * rand(size(R, 1))

rd = RefElemData(Line(), SBP(), N) 
md = MeshData(uniform_mesh(Line(), 8), rd; 
          is_periodic=true)

initial_condition = initial_condition_basic

equations = CompressibleEulerEquations1D(1.4)

u0 = rd.Pq * initial_condition.(md.xq, 0.0, equations)
u = copy(u0)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

Δ = diagm(0 => -ones(rd.N+1), 1=>ones(rd.N+1))[1:end-1,:]
R = tril(ones(size(Δ')), -1)


cache = (; 
           rd, md, 
           B = Diagonal([-1; zeros(rd.N-1); 1]),
           high_order_operators=(; Q_skew = rd.M * rd.Dr - (rd.M * rd.Dr)'), 
           low_order_operators=(; Q_skew_low = Qr-Qr'), 
           fv_operators = (; Δ, R), 
           blend, knapsack_solver = QuadraticKnapsackSolver{Float64}(),
           physics = (; equations, volume_flux = flux_central, surface_flux = flux_lax_friedrichs),
        )

tspan = (0, 2.1)
ode = ODEProblem(rhs!, u, tspan, cache)
integrator = init(ode, SSPRK43()) #; dt = 1e-10)
sol = solve(ode, SSPRK43(), 
            dt = 1 / (md.num_elements * rd.N^2), adaptive=false,
            abstol=1e-10, reltol=1e-8, saveat=LinRange(tspan..., 10); 
            callback=AliveCallback(alive_interval=1000))

rq, wq = gauss_quad(0, 0, N+2)
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x

L1_error = sum(wJq .* map(x -> sum(abs.(x)), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end]))
L2_error = sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end])))
Linf_error = maximum(map(x -> maximum(abs.(x)), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end]))

using Printf
println("N = $N, K1D = $(md.num_elements), blend = $(cache.blend): \n" * 
        "L1_error   = $(@sprintf("%.7f", L1_error)), \n" * 
        "L2_error   = $(@sprintf("%.7f", L2_error)), \n" * 
        "Linf_error = $(@sprintf("%.7f", Linf_error))")

u = sol.u[end]
default()
plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false)
plot!(rd.Vp * md.x, getindex.(initial_condition.(rd.Vp * md.x, tspan[2], equations), 1), leg=false)
