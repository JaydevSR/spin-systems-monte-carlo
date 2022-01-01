include("xywolff.jl")

#=
Perform simulation
=#

N = 20

Temps = [i for i = 0.1:0.2:2.9]
esteps = 5000  # Number of steps for equilibration
nsteps = 10000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

# c_T = zeros(Float64, length(Temps))  # Array of specific heat
# err_c_T = zeros(Float64, length(Temps))

for i = 1:length(Temps)
    global spins
    T = Temps[i]

    E_arr = simulate_xy_wolff(N, T, esteps, nsteps)

    u_T[i] = mean(E_arr) / N^2
    err_u_T[i] = blocking_err(E_arr, A -> mean(A) / N^2)

    # c_T[i] = specific_heat(E_arr, T, N)
    # err_c_T[i] = blocking_err(E_arr, specific_heat, T, N)
end


#=
Plots
=#
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "XY Model for Lattice Size $(N)")

# ax2 = Axis(f[1, 2], xlabel = "temperature, T", ylabel = "specific heat, c",
#     title = "XY Model for Lattice Size $(N)")

errorbars!(
    ax1, Temps, u_T, err_u_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, u_T,
    markersize = 10
)
# errorbars!(
#     ax2, Temps, c_T, err_c_T,
#     whiskerwidth = 10)
# scatter!(ax2, Temps, c_T)

save("XY/XYWolff/plots/internal_energy_vs_Temp_$(N).png", f)

# spins = rand(Float64, (N, N))
# T = 2.0

# for i=1:6000
#     xy_wolff_step!(spins, N, T)
#     if i%500 == 0
#         p = plot_spins(spins, N)
#         display(p)
#     end
# end