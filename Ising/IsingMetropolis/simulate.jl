include("isingmetro.jl")

#=
Perform simulation
=#

N = 20  # Lattice size

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

    E = total_energy(spins)
    M = total_magnetization(spins)

    # Let the system reach equilibrium
    for step = 1:eqsteps
        E, M = ising_metropolis_sweep!(spins, T, E, M)
    end

    u_arr = zeros(Float64, nsteps)
    m_arr = zeros(Float64, nsteps)

    # Iterate for calculating averages
    for step = 1:nsteps
        E, M = ising_metropolis_sweep!(spins, T, E, M)
        u_arr[step] = E / N^2
        m_arr[step] = M / N^2
    end

    u_T[i] = mean(u_arr)
    err_u_T[i] = blocking_err(u_arr, mean)

    m_T[i] = mean(m_arr)
    err_m_T[i] = blocking_err(m_arr, mean)

    c_T[i] = specific_heat(u_arr, T, N)
    err_c_T[i] = blocking_err(u_arr, specific_heat, T, N)

    χ_T[i] = succeptibility(m_arr, T, N)
    err_χ_T[i] = blocking_err(m_arr, succeptibility, T, N)
end


#=
Plots
=#

scatter(Temps, u_T, yerr = u = err_u_T)
xlabel!("temperature, T")
ylabel!("internal energy, u")
title!("Ising Model for Lattice Size $(N)")
savefig("Ising/IsingMetropolis/plots/u_vs_T_$(N).png")

scatter(Temps, m_T, yerr = err_m_T)
xlabel!("temperature, T")
ylabel!("magnetization, m")
title!("Ising Model for Lattice Size $(N)")
savefig("Ising/IsingMetropolis/plots/m_vs_T_$(N).png")

scatter(Temps, c_T, yerr = err_c_T)
xlabel!("temperature, T")
ylabel!("specific heat, c")
title!("Ising Model for Lattice Size $(N)")
savefig("Ising/IsingMetropolis/plots/c_vs_T_$(N).png")

scatter(Temps, χ_T, yerr = err_χ_T)
xlabel!("temperature, T")
ylabel!("succeptibility, χ")
title!("Ising Model for Lattice Size $(N)")
savefig("Ising/IsingMetropolis/plots/x_vs_T_$(N).png")


# Equilibration for some T

# T = 2.0
# energies = Float64[total_energy(spins)]
# magnetizations = Float64[total_magnetization(spins)]

# for step=1:nsteps
#     E, M = ising_metropolis_sweep!(spins, T, energies[step], magnetizations[step])
#     append!(energies, E)
#     append!(magnetizations, M)
# end

# p = plot(1:nsteps, energies[2:end] ./ N^2)
# p = plot!(1:nsteps, magnetizations[2:end] ./ N^2)