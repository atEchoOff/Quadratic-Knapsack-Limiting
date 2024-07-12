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
    (; physics, rd, md, shock_capturing, VDM_inv) = cache
    (; Q_skew) = cache.high_order_operators
    (; Q_skew_low) = cache.low_order_operators    
    (; Δ, R) = cache.fv_operators
    (; equations, volume_flux, surface_flux) = physics

    # interface flux
    uf = rd.Vf * u
    uP = uf[md.mapP]

    interface_flux = @. surface_flux(uf, uP, SVector{1}(md.nxJ), equations)
    du .= rd.LIFT * interface_flux
   
    rhs_vol_high = zeros(eltype(du), size(du, 1), md.num_elements)
    rhs_vol_low = zeros(eltype(du), size(du, 1), md.num_elements)
    
    A = fill!(similar(du, size(du, 1), size(du, 1)), zero(eltype(du)))
    a = zeros(Float64, size(du, 1) + 1, md.num_elements)
    b = zero(Float64)

    for e in 1:md.num_elements
        rhs_vol_high_element = view(rhs_vol_high, :, e)
        rhs_vol_low_element = view(rhs_vol_low, :, e)
        
        # high order update
        u_element = view(u, :, e)
        a_element = view(a, :, e)

        v = cons2entropy.(u_element, equations)
        for i in eachindex(u_element), j in eachindex(u_element)            
            if i > j
                fij = volume_flux(u_element[i], u_element[j], 1, equations)

                FH_ij = Q_skew[i, j] * fij
                rhs_vol_high_element[i] +=  FH_ij
                rhs_vol_high_element[j] += -FH_ij

                FL_ij = zero(eltype(rhs_vol_low_element))
                if abs(Q_skew_low[i, j]) > 1e-12
                    nij = Q_skew_low[i, j] / abs(Q_skew_low[i, j])
                    fij = surface_flux(u_element[i], u_element[j], SVector{1}(nij), equations)
                    FL_ij = abs(Q_skew_low[i, j]) * fij
                    rhs_vol_low_element[i] += FL_ij
                    rhs_vol_low_element[j] -= FL_ij    
                end

                # create antidiffusive flux matrix
                A[i,j] = FH_ij - FL_ij
                A[j,i] = -A[i, j]
            end
        end
        
        # optimization target
        my_b = sum(cache.B * psi.(u_element, equations)) + sum(dot.(v, rhs_vol_low_element))

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

            my_b = (1 - .01 * ϵ) * my_b
        end

        b += my_b

        a_element .= (v' * Δ * Diagonal(R * (rhs_vol_low_element - rhs_vol_high_element)))'
    end

    θ = cache.knapsack_solver(vec(a), b)
    θ = reshape(θ, (size(du, 1) + 1, md.num_elements))

    # Set du for each element
    for e in 1:md.num_elements
        rhs_vol_high_element = view(rhs_vol_high, :, e)
        rhs_vol_low_element = view(rhs_vol_low, :, e)
        θ_element = view(θ, :, e)
        view(du, :, e) .+= rd.M \ (Δ * (Diagonal(θ_element) * R * rhs_vol_high_element + Diagonal(1 .- θ_element) * R * rhs_vol_low_element))
    end
    
    # invert Jacobian and mass matrix
    du ./= -md.J
end

psi(u, ::CompressibleEulerEquations1D) = u[2] # rho * v1

function initial_condition_basic(x, t, equations::CompressibleEulerEquations1D)           
    # rho = 1.0 + exp(-100 * (x-0.25)^2)
    # # rho = 1.0 + .75 * (x > 0)
    # u = 0.0
    # p = rho^equations.gamma 

    u = 1.0
    rho = 1 + .999 * sin(pi * (x - u * t))
    p = 0.1

    return SVector(prim2cons(SVector(rho, u, p), equations))
end

initial_condition = initial_condition_basic

initial_limiting_coefficients(R, u; delta=0) = 1 .- delta * rand(size(R, 1))

rd = RefElemData(Line(), SBP(), N) 
(VX, ), EToV = uniform_mesh(Line(), K)
md = MeshData((VX,), EToV, rd)
md = make_periodic(md)

equations = CompressibleEulerEquations1D(1.4)

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
           blend = :subcell,
           knapsack_solver = knapsack_solver,
        physics = (; equations, volume_flux = volume_flux, surface_flux = flux_lax_friedrichs),
        shock_capturing = shock_capturing,
        VDM_inv = inv(rd.VDM)
        )

# dt = 0.25 * estimate_h(rd, md) * (1 / (2 * rd.N + 1))
tspan = (0, 2.5)
ode = ODEProblem(rhs!, u, tspan, cache)
@time sol = solve(ode, timestepper, dt = dt, abstol=abstol, reltol=reltol, saveat=LinRange(tspan..., 10), adaptive=adaptive)

# @gif for u in sol.u
#     plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false, ylims=(0, 2))
# end
u = sol.u[end]
plot(rd.Vp * md.x, rd.Vp * getindex.(u, 1), leg=false)
# plot!(rd.Vp * md.x, getindex.(initial_condition.(rd.Vp * md.x, tspan[2], equations), 1), leg=false)

rq, wq = gauss_quad(0, 0, N+2)
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM
wJq = Diagonal(wq) * (Vq * md.J)
xq = Vq * md.x

L1_error = sum(wJq .* map(x -> sum(abs.(x)), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end]))
L2_error = sqrt(sum(wJq .* map(x -> sum(x.^2), initial_condition.(xq, sol.t[end], equations) - Vq * sol.u[end])))

println("N = $N, K1D = $(md.num_elements), L1_error = $L1_error, L2_error = $L2_error")