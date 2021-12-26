#=
Implementaion of the Wolff algorithm for Monte-Carlo simulation of the Ising model.
The Hamiltonian has B=0 and J=1. Also, k=1 => β = 1/T. 
=#
using Plots, Statistics

function ising_wolff_step!(spins, P_add)
    N = size(spins)[1]
    seed = rand(1:N, 2)  # seed spin position
    cluster = get_cluster!(spins, seed, P_add)

    ΔM = -2*spins[seed...]*sum(cluster)
    # flip cluster
    @. spins = ifelse(cluster, -spins, spins)
    return ΔM
end


function get_cluster!(spins, seed, P_add)
    cluster = falses(size(spins))
    sval = spins[seed...]
    stack = [seed]
    cluster[seed...] = true
    while !isempty(stack)
        k = pop!(stack)
        for δ ∈ ([1, 0], [N-1, 0], [0, 1], [0, N-1])
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            if spins[nn...] == sval && !cluster[nn...] && rand() < P_add
                push!(stack, nn)
                cluster[nn...] = true
            end
        end
    end
    return cluster
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


"""
    autocorrelation_fn(mags, N)

Calculate the autocorrelation function (normalized) of the given time series array.
"""
function autocorrelation_fn(series, N)
    tmax = length(series)
    autocorr = zeros(Float64, tmax)
    for t ∈ 1:tmax-1
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for tk ∈ 1:tmax-t
            sum1 += series[tk]*series[tk+t]
            sum2 += series[tk]
            sum3 += series[tk+t]
        end
        autocorr[t] = sum1 / (tmax-t) - (sum2*sum3) / (tmax - t)^2
    end
    @. autocorr /= N^2
    @. autocorr /= autocorr[1]
    return autocorr
end


"""
    bootstrap_err(samples, calc_qty; r=100)

Estimate the error in the given samples by bootstrap method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `r` is a keyword arguments giving number of resamples.
"""
function bootstrap_err(samples, calc_qty, args...; r=100)
    nob = length(samples)
    resample_arr = zeros(Float64, nob)
    for i=1:r
        resample = rand(samples, nob)
        resample_arr[i] = calc_qty(resample, args...)
    end
    err = std(resample_arr, corrected=false)
    return err
end


"""
    blocking_err(samples, calc_qty; blocks=10)

Estimate the error in the given samples by blocking method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `blocks` is a keyword arguments giving number of blocks.
"""
function blocking_err(samples, calc_qty, args...; blocks=20)
    block_array = zeros(Float64, blocks)
    blocklength = length(samples) ÷ blocks
    for i=1:blocks
        sample_block = samples[(i-1)*blocklength + 1 : i*blocklength]
        block_array[i] = calc_qty(sample_block, args...)
    end
    err = std(block_array)
    return err
end


"""
    specific_heat(u_vals, T, N)

Calculate the specific heat from given array of internal energy per site (`N²` sites) at temperature `T`.
"""
function specific_heat(u_vals, T, N)
    return (T^-2) * N^2 * var(u_vals, corrected=false)
end


"""
    succeptibility(m_vals, T, N)

Calculate the succeptibility from given array of mean magnetization per site (`N²` sites) at temperature `T`.
"""
function succeptibility(m_vals, T, N)
    return (T^-2) * N^2 * var(m_vals, corrected=false)
end


#=
Perform simulation
=#

N = 10  # Lattice size

# Initialize lattice
spins = ones(N, N)  # T = 0
# spins = rand([-1, 1], (N, N))  # T = ∞


Temps = [i for i=1.0:0.1:3.6]
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

for i=1:length(Temps)
    global spins
    T = Temps[i]
    P_add = 1 - exp(-2/T)

    # Let the system reach equilibrium
    for step=1:eqsteps
        ΔM = ising_wolff_step!(spins, P_add)
    end

    u_arr = zeros(Float64, nsteps)
    m_arr = zeros(Float64, nsteps)

    u_arr[1] = total_energy(spins) / N^2
    m_arr[1] = total_magnetization(spins) / N^2


    # Iterate for calculating averages
    for step=1:nsteps-1
        ΔM = ising_wolff_step!(spins, P_add)
        u_arr[step+1] = total_energy(spins) / N^2
        m_arr[step+1] = ΔM / N^2
    end

    m_arr = abs.(m_arr)
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

scatter(Temps, u_T, yerr=u=err_u_T)
xlabel!("temperature, T")
ylabel!("internal energy, u")
title!("Ising Model for Lattice Size $(N)")
savefig("IsingWolff/plots/u_vs_T_$(N).png")

scatter(Temps, m_T, yerr=err_m_T)
xlabel!("temperature, T")
ylabel!("magnetization, m")
title!("Ising Model for Lattice Size $(N)")
savefig("IsingWolff/plots/m_vs_T_$(N).png")

scatter(Temps, c_T, yerr=err_c_T)
xlabel!("temperature, T")
ylabel!("specific heat, c")
title!("Ising Model for Lattice Size $(N)")
savefig("IsingWolff/plots/c_vs_T_$(N).png")

scatter(Temps, χ_T, yerr=err_χ_T)
xlabel!("temperature, T")
ylabel!("succeptibility, χ")
title!("Ising Model for Lattice Size $(N)")
savefig("IsingWolff/plots/x_vs_T_$(N).png")


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