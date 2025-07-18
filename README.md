# Quadratic-Knapsack-Limiting

This is the testing repo for enforcing a cell entropy inequality with quadratic knapsack limiting. In order to run tests, first run either 1D/Running_Interface.jl or 2D/Running_Interface.jl. Here is an explanation of each tunable setting:

- `N` refers to the degree of polynomial approximation, and it should be a positive integer.
- `K` refers to the number of elements along each individual dimension. It should be a positive integer. 
- `knapsack_solver` refers to the chosen algorithm for knapsack limiting. Options include:
    - `QuadraticKnapsackMinimizer{Float64}()` performs quadratic knapsack limiting.
    - `ContinuousKnapsackSolver((N + 1) * (N + 2))` performs linear knapsack limiting.
    - `QuadraticKnapsackSolverA{Float64}()` performs weighted quadratic knapsack limiting.
    - `NonKnapsackSolver{Float64}()` sets all blending coefficients to 0.
- `volume_flux` refers to the chosen volume flux in the scheme. By default, users should choose `flux_central` for knapsack limiting. However, to test entropy conservative flux differencing, you may choose `NonKnapsackSolver` as the blending strategy, and a valid entropy conservative `volume_flux`. 
- `blend` refers to the blending strategy type. Note that quadratic knapsack limiting is a *minimization* problem, whereas linear knapsack limiting was originally framed as a *maximization* problem. Valid options include:
    - `subcell` for linear knapsack limiting, or for entropy conservative flux differencing in combination with `NonKnapsackSolver`.
    - `subcell_reversed` for quadratic knapsack limiting
    - `viscosity` for testing a discrete artificial viscosity scheme (no knapsack)
    - `elementwise` for elementwise-constant blending coefficients (no knapsack)
    - `loworder` for the purely low order method (no knapsack)
    - Some additional strategies exist, such as `subcell_dynamic` which applies subcell with different blending coefficients for different conserved variables, but these are not tested and not used in research.
- `shock_capturing` applies a basic shock capturing to `subcell` schemes. It is not used in research, and setting it to 0 disables it. It should be a value between 0 and 1.
- `nodewise_shock_capturing` also applies a subcell shock capturing in `subcell` schemes, but it is also not used in research, and setting it to 0 disables it. 
- `abstol` and `reltol` define absolute and relative tolerances for adaptive timestepping.
- `timestepper` is the chosen timestepper.
- `adaptive` defines whether or not to use adaptive timestepping
- `dt` defines the fixed timestep for `adaptive = false`, or the initial timestep for `adaptive = true`. 
- `saveat` defines the timesteps at which to save the resulting solution.
- `preserve_positivity` refers to the chosen relative positivity condition. Setting it to -1 disables it, setting all limiting coefficients to 0 (allowing blending coefficients to vary between 0 and 1).

After running 1D(or 2D)/RunningInterface.jl, you may run whichever problem you want. There are additional files, such as `ErrorPlotterScript`, `CFL_over_time`, `AdaptiveTimesteppingScript`, or `SpaceOrderScript` available for analyzing the resulting data. 