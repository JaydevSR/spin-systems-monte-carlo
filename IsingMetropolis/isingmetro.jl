#=
Implementaion of the Metropolis algorithm for Monte-Carlo simulation of the Ising model.
The Hamiltonian has B=0 and J=1. Also, k=1 => β = 1/T. 
=#

using Plots


"""
    ising_metropolis_sweep(spins, T, E, M)

Perform one sweep of the lattice using single-spin-flip dynamics (1 sweep == N*N flip attempts).
Here arguments E and M are the total energy and total magnetization before the sweep.
Returns total energy and magnetization after sweep.
"""
function ising_metropolis_sweep!(spins, T, E, M)
    N = size(spins)[1]
    for i=1:N^2
        k = rand(1:N, 2)
        ΔE = deltaE(spins, k, N)
        if accept_flip(ΔE, T)
            spins[k[1], k[2]] *= -1
            E = E + ΔE
            M = M + 2spins[k[1], k[2]]
        end
    end
    return E, M
end


"""
    accept_flip(spins, N, k_i, k_j, T)

Determine whether to accept the next state or not according to Metropolis acceptance rates.
Returns `true` or `false`.
"""
function accept_flip(ΔE, T)
    # Metropolis Acceptance Rates
    if ΔE <= 0
        return true
    elseif rand() < exp(-ΔE/T)
        return true
    else
        return false
    end
end


"""
    delta_E(spins, k)

Calculate the energy difference between two states for one spin flip at site k.
"""
function deltaE(spins, k, N)
    ΔE = 0
    for δ ∈ ([1, 0], [N-1, 0], [0, 1], [0, N-1])
        nn = k + δ
        @. nn = mod1(nn, N)  # Apply periodic boundary conditions
        ΔE += spins[nn[1], nn[2]]
    end
    ΔE *= 2spins[k[1], k[2]]
end


"""
    total_magnetization(spins)

Calculate the total magnetization of square spin lattice.
"""
function total_magnetization(spins)
    return sum(sum(spins))
end


"""
    total_energy(spins)

Calculate the total energy of the square spin lattice (with zero field and J=1).
"""
function total_energy(spins)
    N = size(spins)[1]
    running_sum = 0
    for i=1:N
        for j=1:N
            s_k = spins[i, j]
            for δ ∈ ([1, 0], [N-1, 0], [0, 1], [0, N-1])
                nn = [i, j] + δ
                @. nn = mod1(nn, N)  # Apply periodic boundary conditions
                running_sum += s_k*spins[nn[1], nn[2]]
            end
        end
    end
    return -running_sum / 2  # divide by 2 because each bond counted twice
end


#=
Perform simulation
=#

N = 10  # Lattice size

# Initialize lattice
# spins = ones(N, N)  # T = 0
spins = rand([-1, 1], (N, N))  # T = ∞


Temps = [i for i=1.0:0.1:3.6]
eqsteps = 5000  # Number of steps for equilibration
nsteps = 10000  # Number of steps for measurements

u_T = zeros(Float32, length(Temps))  # Array of internal energy per site
m_T = zeros(Float32, length(Temps))  # Array of mean magnetization per site
c_T = zeros(Float32, length(Temps))  # Array of specific heat
χ_T = zeros(Float32, length(Temps))  # Array of succeptibility

for i=1:length(Temps)
    global spins
    T = Temps[i]
    E = total_energy(spins)
    M = total_magnetization(spins)
    E_sum = M_sum = E_sq_sum = M_sq_sum = 0

    # Let the system reach equilibrium
    for step=1:eqsteps
        E, M = ising_metropolis_sweep!(spins, T, E, M)
    end

    # Iterate for calculating averages
    for step=1:nsteps
        E, M = ising_metropolis_sweep!(spins, T, E, M)
        E_sum += E
        M_sum += M
        E_sq_sum += E*E
        M_sq_sum += M*M
    end

    u_T[i] = ( E_sum / (nsteps*N*N) )
    m_T[i] = ( M_sum / (nsteps*N*N) )
    c_T[i] = (T^-2) * ( E_sq_sum / (nsteps*N*N) - E_sum*E_sum / (nsteps*nsteps*N*N) )
    χ_T[i] = (T^-1) * ( M_sq_sum / (nsteps*N*N) - M_sum*M_sum / (nsteps*nsteps*N*N) )
end


#=
Plots
=#

scatter(Temps, u_T)
xlabel!("temperature, T")
ylabel!("internal energy, u")
savefig("IsingMetropolis/plots/u_vs_T.png")

scatter(Temps, m_T)
xlabel!("temperature, T")
ylabel!("magnetization, m")
savefig("IsingMetropolis/plots/m_vs_T.png")

scatter(Temps, c_T)
xlabel!("temperature, T")
ylabel!("specific heat, c")
savefig("IsingMetropolis/plots/c_vs_T.png")

scatter(Temps, χ_T)
xlabel!("temperature, T")
ylabel!("succeptibility, χ")
savefig("IsingMetropolis/plots/x_vs_T.png")


# T = 2.0
# energies = Float32[total_energy(spins)]
# magnetizations = Float32[total_magnetization(spins)]

# for step=1:nsteps
#     E, M = ising_metropolis_sweep!(spins, T, energies[step], magnetizations[step])
#     append!(energies, E)
#     append!(magnetizations, M)
# end

# p = plot(1:nsteps, energies[2:end] ./ N^2)
# p = plot!(1:nsteps, magnetizations[2:end] ./ N^2)