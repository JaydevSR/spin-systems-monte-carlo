using CairoMakie

include("xywolff.jl")

#=
Perform simulation
=#

N = 10
println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

Temps = [i for i = 0.1:0.2:2.9]
esteps = 5000  # Number of steps for equilibration
nsteps = 50000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

c_T = zeros(Float64, length(Temps))  # Array of specific heat
err_c_T = zeros(Float64, length(Temps))

for i = 1:length(Temps)
    global spins
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    E_arr = simulate_xy_wolff(N, T, esteps, nsteps)

    u_T[i] = mean(E_arr) / N^2
    err_u_T[i] = blocking_err(E_arr, A -> mean(A) / N^2)

    c_T[i] = specific_heat(E_arr, T, N)
    err_c_T[i] = blocking_err(E_arr, specific_heat, T, N)
    println("   |          ")
    println("   +-> Done.\n")
end


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "XY Model for Lattice Size $(N)")

ax2 = Axis(f[1, 2], xlabel = "temperature, T", ylabel = "specific heat, c",
    title = "XY Model for Lattice Size $(N)")

errorbars!(
    ax1, Temps, u_T, err_u_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, u_T,
    markersize = 10
)

errorbars!(
    ax2, Temps, c_T, err_c_T,
    whiskerwidth = 10)
scatter!(
    ax2, Temps, c_T,
    markersize=10
)

save("XY/XYWolff/plots/u&c_vs_T_$(N).png", f)

# spins = rand(Float64, (N, N))
# T = 2.0

# for i=1:6000
#     xy_wolff_step!(spins, N, T)
#     if i%500 == 0
#         p = plot_spins(spins, N)
#         display(p)
#     end
# end
println("Program Finished!")
println("===========================\n")