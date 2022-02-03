include("../src/spinmc.jl")

#=
Perform simulation
=#

N = 20
println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

Temps = append!([i for i = 0.1:0.1:0.7], [i for i = 0.8:0.05:1.4], [i for i = 1.5:0.1:2.5])
esteps = 2000  # Number of steps for equilibration
nsteps = 20000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

c_T = zeros(Float64, length(Temps))  # Array of specific heat
err_c_T = zeros(Float64, length(Temps))

for i = 1:length(Temps)
    global spins
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    u_arr = simulate_xy_wolff(N, T, esteps, nsteps, from_infinity = false) ./ N^2

    u_T[i] = mean(u_arr) 
    err_u_T[i] = blocking_err(u_arr, A -> mean(A))

    c_T[i] = specific_heat(u_arr, T, N)
    err_c_T[i] = blocking_err(u_arr, specific_heat, T, N)
    println("   |          ")
    println("   +-> Done.\n")
end


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "XY $(N)x$(N): Internal energy (per site) v/s temperature")

ax2 = Axis(f[2, 1], xlabel = "temperature, T", ylabel = "specific heat, c",
    title = "XY $(N)x$(N): Specific heat (per site) v/s temperature")

errorbars!(
    ax1, Temps, u_T, err_u_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, u_T,
    markersize = 7
)

errorbars!(
    ax2, Temps, c_T, err_c_T,
    whiskerwidth = 10)
scatter!(
    ax2, Temps, c_T,
    markersize = 7
)

save("simulations/results/xy/u_and_c_vs_Temp_$(N)_wolff.png", f)

println("Program Finished!")
println("===========================\n")