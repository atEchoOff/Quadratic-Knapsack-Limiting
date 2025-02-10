include("ObjectTransitioner.jl")

println("This terminal is running case $(case)")

testcases = (
    (; N = 3, K = 64, timestepper = RK4(), initial_condition = initial_condition_sod_shock),
    (; N = 7, K = 32, timestepper = RK4(), initial_condition = initial_condition_sod_shock),
    (; N = 3, K = 64, timestepper = SSPRK43(), initial_condition = initial_condition_sod_shock),
    (; N = 7, K = 32, timestepper = SSPRK43(), initial_condition = initial_condition_sod_shock),
    (; N = 3, K = 64, timestepper = Tsit5(), initial_condition = initial_condition_sod_shock),
    (; N = 7, K = 32, timestepper = Tsit5(), initial_condition = initial_condition_sod_shock)
)

for testcase in testcases
    global (; N, K, timestepper, initial_condition) = testcase

    include("1D/ErrorPlotterScript.jl")

    save("N$(N)K$(K)TS$(nameof(typeof(timestepper)))C$(case)", errors)
end