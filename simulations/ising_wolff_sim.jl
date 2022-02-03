include("../src/spinmc.jl")

#=
Perform simulation
=#

N = 20  # Lattice size
println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

# Initialize lattice
spins = ones(N, N)  # T = 0
# spins = rand([-1, 1], (N, N))  # T = ∞


Temps = [i for i = 1.0:0.1:3.6]
eqsteps = 2000  # Number of steps for equilibration
nsteps = 6000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

m_T = zeros(Float64, length(Temps))  # Array of mean magnetization per site
err_m_T = zeros(Float64, length(Temps))

c_T = zeros(Float64, length(Temps))  # Array of specific heat
err_c_T = zeros(Float64, length(Temps))

χ_T = zeros(Float64, length(Temps))  # Array of succeptibility
err_χ_T = zeros(Float64, length(Temps))

for i = 1:length(Temps)
    global spins
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    P_add = 1 - exp(-2 / T)

    # Let the system reach equilibrium
    for step = 1:eqsteps
        ΔM = isingwolff_step!(spins, P_add)
    end

    u_arr = zeros(Float64, nsteps)
    m_arr = zeros(Float64, nsteps)

    u_arr[1] = ising_total_energy(spins) / N^2
    m_arr[1] = ising_total_magnetization(spins) / N^2


    # Iterate for calculating averages
    for step = 1:nsteps-1
        ΔM = isingwolff_step!(spins, P_add)
        u_arr[step+1] = ising_total_energy(spins) / N^2
        m_arr[step+1] = ising_total_magnetization(spins) / N^2
    end

    m_arr = abs.(m_arr)  # take absolute magnetization
    u_T[i] = mean(u_arr)
    err_u_T[i] = blocking_err(u_arr, mean)

    m_T[i] = mean(m_arr)
    err_m_T[i] = blocking_err(m_arr, mean)

    c_T[i] = specific_heat(u_arr, T, N)
    err_c_T[i] = blocking_err(u_arr, specific_heat, T, N)

    χ_T[i] = succeptibility(m_arr, T, N)
    err_χ_T[i] = blocking_err(m_arr, succeptibility, T, N)
    println("   |          ")
    println("   +-> Done.\n")
end


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(
    f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "Ising $(N)x$(N): Internal energy vs temperature"
    )

ax2 = Axis(
    f[1, 2], xlabel = "temperature, T", ylabel = "specific heat, c",
    title = "Ising $(N)x$(N): Specific heat vs temperature"
    )

ax3 = Axis(
    f[2, 1], xlabel = "temperature, T", ylabel = "magnetization, m",
    title = "Ising $(N)x$(N): Magnetization vs temperature"
    )

ax4 = Axis(
    f[2, 2], xlabel = "temperature, T", ylabel = "succeptibility, χ",
    title = "Ising $(N)x$(N): Succeptibility vs temperature"
    )


errorbars!(
    ax1, Temps, u_T, err_u_T,
    whiskerwidth = 10
)
scatter!(ax1, Temps, u_T)


errorbars!(
    ax2, Temps, c_T, err_c_T,
    whiskerwidth = 10
)
scatter!(ax2, Temps, c_T)

errorbars!(
    ax3, Temps, m_T, err_m_T,
    whiskerwidth = 10
)
scatter!(ax3, Temps, m_T)


errorbars!(
    ax4, Temps, χ_T, err_χ_T,
    whiskerwidth = 10
)
scatter!(ax4, Temps, χ_T)

save("simulations/results/ising/physical_qtys_vs_Temp_$(N)_wolff.png", f)

println("Program Finished!")
println("===========================\n")