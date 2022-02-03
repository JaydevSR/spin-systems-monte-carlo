"""
    ising_total_magnetization(spins)

Calculate the total magnetization of square spin lattice.
"""
function ising_total_magnetization(spins)
    return sum(spins)
end


"""
    ising_total_energy(spins)

Calculate the total energy of the square spin lattice (with zero field and J=1).
"""
function ising_total_energy(spins)
    N = size(spins)[1]
    running_sum = 0
    for i = 1:N
        for j = 1:N
            s_k = spins[i, j]
            for δ ∈ ([1, 0], [N - 1, 0], [0, 1], [0, N - 1])
                nn = [i, j] + δ
                @. nn = mod1(nn, N)  # Apply periodic boundary conditions
                running_sum += s_k * spins[nn[1], nn[2]]
            end
        end
    end
    return -running_sum / 2  # divide by 2 because each bond counted twice
end

"""
    xy_total_energy(spins)

Calculate the total energy of the square spin lattice `spins` of size `(N, N)` with zero field and J=1.
"""
function xy_total_energy(spins::Matrix, N::Int64)
    running_sum = 0
    for idx in CartesianIndices(spins)
        s_k = spins[idx]
        for δ ∈ CartesianIndex.([(1, 0), (N - 1, 0), (0, 1), (0, N - 1)])
            nn = idx + δ
            nn = CartesianIndex(mod1.(Tuple(nn), N))  # Apply periodic boundary conditions
            running_sum += cos2pi(s_k - spins[nn])
        end
    end
    return - running_sum / 2  # divide by 2 because each bond counted twice
end

"""
    specific_heat(u_vals, T, N)

Calculate the specific heat from given array of internal energy per site (`N²` sites) at temperature `T`.
"""
function specific_heat(u_vals, T, N)
    return (T^-2) * N^2 * var(u_vals, corrected = false)
end


"""
    succeptibility(m_vals, T, N)

Calculate the succeptibility from given array of mean magnetization per site (`N²` sites) at temperature `T`.
"""
function succeptibility(m_vals, T, N)
    return (T^-2) * N^2 * var(m_vals, corrected = false)
end