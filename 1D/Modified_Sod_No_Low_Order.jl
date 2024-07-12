using StartUpDG
using StaticArrays, StructArrays, RecursiveArrayTools
using LinearAlgebra
using Trixi
using OrdinaryDiffEq
using Plots
using CSV
using Tables

include("../L1_knapsack.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_weighted_a.jl")

function rhs!(du, u, cache, t)
    (; physics, rd, md) = cache
    (; Q_skew) = cache.high_order_operators
    (; Q_skew_low) = cache.low_order_operators    
    (; Δ, R) = cache.fv_operators
    (; equations, volume_flux, surface_flux, gamma) = physics
    (; bc, VDM_inv, shock_capturing, nodewise_shock_capturing) = cache

    # Boundary conditions (set end points to initial conditions)
    u[1, 1] = bc[:, 1]
    u[end, end] = bc[:, 2]

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
                ## Chan's flux function
                fij = volume_flux(u_element[i], u_element[j], 1, equations)
                ##
                ## Raymond's flux function
                # fij = volume_flux(u_element[i], u_element[j], gamma);
                ##
                FH_ij = Q_skew[i, j] * fij
                rhs_vol_high[i] +=  FH_ij
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
        
        # optimization target
        b = sum(cache.B * psi.(u_element, equations))

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

            b = (1 - .02 * ϵ) * b
        end

        (; blend) = cache

        if blend == :elementwise
            # 1. elementwise blending of rhs_vol_high/low
            denom = (sum(dot.(v, rhs_vol_high - rhs_vol_low)) + 100 * eps())
            θ = b / denom
            θ = max(0.0, min(1.0, θ))
            view(du, :, e) .+= rd.M \ (θ * rhs_vol_high + (1 - θ) * rhs_vol_low)
        elseif blend == :subcell
            a = v' * Δ * Diagonal(R * rhs_vol_high)


            # Call the Knapsack Solver
            a = vec(a)
            θ = cache.knapsack_solver(a, b)
            
            if nodewise_shock_capturing
                epsilon = 240.0
                # @. a *= -log(1/e * θ)
                # @. a *= 1 / cbrt(θ)
                @. a *= 1 - 4 * N^2 * log(min(θ + 1 / (2N * K), 1))
                # # @. a *= (1 - 2 * N * (K / 64) * log(θ + 1e-8))
                # @. a *= -epsilon * θ + (1 + epsilon)
                # # @. a *= (-sqrt(epsilon) * θ + sqrt(epsilon))^2 + 1
                # @. a *= -epsilon * min(1, θ + 1e-3) + (1 + epsilon)

                θ = cache.knapsack_solver(a, b)
                # @. a *= -epsilon * θ + (1 + epsilon)
                # # @. a *= (-sqrt(epsilon) * θ + sqrt(epsilon))^2 + 1

                # θ = cache.knapsack_solver(a, b)
            end
            view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ) * R * rhs_vol_high))
        end
    end
    
    # invert Jacobian and mass matrix
    du ./= -md.J
end

psi(u, ::CompressibleEulerEquations1D) = u[2] # rho * v1

function initial_condition_modified_sod(x, t, equations::CompressibleEulerEquations1D)
    if x[1] < .2
        rho = 1.0
        v1 = .75
        p = 1.0
    else
        rho = .125
        v1 = 0.0
        p = .1
    end
    return prim2cons(SVector(rho, v1, p), equations)
end

function domain_change(x)
    a = 0.0;
    b = 1.0;
    return (b - a)/2 * (x + 1) + a;
end

### Set Initial conditions
initial_condition = initial_condition_modified_sod

rd = RefElemData(Line(), SBP(), N) 
(VX, ), EToV = uniform_mesh(Line(), K)

md = MeshData((VX,), EToV, rd)

gamma = 1.4
equations = CompressibleEulerEquations1D(gamma)

x = domain_change.(rd.Vp * md.x)
x0 = domain_change.(md.xq)
u0 = rd.Pq * initial_condition.(x0, 0.0, equations)
u = copy(u0)

# low order operators
(Qr,), E = sparse_low_order_SBP_operators(rd)

# subcell FV operators
Δ = diagm(0 => -ones(rd.N+1), 1=>ones(rd.N+1))[1:end-1,:]
R = tril(ones(size(Δ')), -1)

# Build the cache for RHS!
cache = (; 
           rd, md, 
           B = Diagonal([-1; zeros(rd.N-1); 1]),
           high_order_operators=(; Q_skew = rd.M * rd.Dr - (rd.M * rd.Dr)'), 
           low_order_operators=(; Q_skew_low = Qr-Qr'), 
           fv_operators = (; Δ, R),
           bc = [u0[1, 1] u0[end, end]],
           VDM_inv = inv(rd.VDM),
           shock_capturing = shock_capturing,
           nodewise_shock_capturing = nodewise_shock_capturing,
           blend = blend,
        knapsack_solver = knapsack_solver,
        physics = (; equations, volume_flux = volume_flux, surface_flux = flux_lax_friedrichs, gamma = gamma),
        );

tspan = (0, .4)
ode = ODEProblem(rhs!, u, tspan, cache)

sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, callback=AliveCallback(alive_interval=1000), adaptive=adaptive)

println("Completed run with N = $N, K = $K, knapsack_solver = $(typeof(knapsack_solver)), timestepper = $(typeof(timestepper)), abstol = $abstol, reltol = $reltol")

u = sol.u[end]
u_plot = rd.Vp * getindex.(u, 1)
plot(x, u_plot, leg=false)

# @gif for u in sol.u
#     plot(x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0, 1))
# end