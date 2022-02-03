include("../src/spinmc.jl")

#=
Perform simulation
=#

N = 10 # [10, 20, 30]  # Lattice size

println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

# Initialize lattice
spins = ones(N, N)  # T = 0

Temps = [i for i = 1.0:0.1:3.6]
eqsteps = 2000  # Number of steps for equilibration
nsteps = 6000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

c_T = zeros(Float64, length(Temps))  # Array of specific heat
err_c_T = zeros(Float64, length(Temps))

for i = 1:length(Temps)
    global spins
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    E = ising_total_energy(spins)
    M = ising_total_magnetization(spins)

    # Let the system reach equilibrium
    for step = 1:eqsteps
        E, M = isingmetro_step!(spins, T, E, M)
    end

    u_arr = zeros(Float64, nsteps)

    # Iterate for calculating averages
    for step = 1:nsteps
        E, M = isingmetro_step!(spins, T, E, M)
        u_arr[step] = E / N^2
    end

    u_T[i] = mean(u_arr)
    err_u_T[i] = blocking_err(u_arr, mean)

    c_T[i] = specific_heat(u_arr, T, N)
    err_c_T[i] = blocking_err(u_arr, specific_heat, T, N)
    println("   |          ")
    println("   +-> Done.\n")
end

Temps2 = [i for i = 1.0:0.01:3.6]
u_T_exact = ising_exact_u.(Temps2)
c_T_exact = ising_exact_c.(Temps2)


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(
    f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "Ising Model for Lattice Size $(N)"
    )

ax2 = Axis(
    f[1, 2], xlabel = "temperature, T", ylabel = "specific heat, c",
    title = "Ising Model for Lattice Size $(N)"
    )

# internal energy per site

lines!(ax1, Temps2, u_T_exact, label = "Exact Solution")

errorbars!(
    ax1, Temps, u_T, err_u_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, u_T,
    markersize = 10, label = "MC results"
)

axislegend(ax1, orientation = :horizontal, position = :lt)

# specific heat per site

lines!(ax2, Temps2, c_T_exact, label = "Exact Solution")

errorbars!(
    ax2, Temps, c_T, err_c_T,
    whiskerwidth = 10)
    
scatter!(
    ax2, Temps, c_T,
    markersize = 10, label = "MC results"
    )

axislegend(ax2, orientation = :horizontal, position = :lt)

println("Saving Plots ...")
save("simulations/results/ising/u_and_c_vs_Temp_$(N)_exact_vs_mc.png", f)

println("Program Finished!")
println("===========================\n")