using Statistics

# TODO : Program not working as exected, fix it!


"""
    xy_wolff_step!(spins::Matrix, N::Int64, T::Float64)

Perform one step of Wolff algorithm for XY model (lattice `spins` of size `(N, N)` at temperature `T`).
"""
function xy_wolff_step!(spins::Matrix, N::Int64, T::Float64)
    u_flip = rand()  # Random unit vector in xy plane
    isingpll, isingprp, lenpll, lenprp = decompose_XY_to_ising(spins, u_flip)
    xy_ising_wolff_step!(isingpll, lenpll, N, T)
    spins = compose_ising_to_XY(isingpll, isingprp, lenpll, lenprp, u_flip)
end

function xy_ising_wolff_step!(isingpll, lenpll, N, T)
    seed = rand(1:N, 2)  # seed spin position
    cluster = get_ising_cluster!(isingpll, lenpll, seed, T)
    # flip cluster
    @. isingpll = ifelse(cluster, -isingpll, isingpll)
end


function get_ising_cluster!(ising_sign, ising_len, seed, T)
    cluster = falses(size(ising_sign))
    sval = ising_sign[seed...]
    stack = [seed]
    cluster[seed...] = true
    while !isempty(stack)
        k = pop!(stack)
        for δ ∈ ([1, 0], [N - 1, 0], [0, 1], [0, N - 1])
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            P_add = 1 - exp(- 2 * ising_len[seed...] * ising_len[nn...] / T)
            if ising_sign[nn...] == sval && !cluster[nn...] && rand() < P_add
                push!(stack, nn)
                cluster[nn...] = true
            end
        end
    end
    return cluster
end

function decompose_XY_to_ising(XYlattice, u_flip)
    along_u = cos2pi.(XYlattice .- u_flip)
    normal_u = sin2pi.(XYlattice .- u_flip)
    Ising_pll = sign.(along_u)
    Ising_prp = sign.(normal_u)
    len_pll = abs.(along_u)
    len_prp = abs.(normal_u)
    return Ising_pll, Ising_prp, len_pll, len_prp
end

function compose_ising_to_XY(ising_pll, ising_prp, len_pll, len_prp, u_flip)
    along_u = ising_pll .* len_pll
    normal_u = ising_prp .* len_prp
    x = along_u .* cos2pi(u_flip) - normal_u .* sin2pi(u_flip)
    y = along_u .* sin2pi(u_flip) + normal_u .* cos2pi(u_flip)
    XYlattice = mod1.(atan.(y, x) ./ (2pi), 1)
    return XYlattice
end


"""
    P_add(u_flip::Float64, s1::Float64, s2::Float64, T::Float64)

Calculate the probability of adding spin `s2` to cluster of as a neighbour of `s1` at temperature `T` w.r.t angle `u_flip`.
"""
function P_add(u_flip::Float64, s1::Float64, s2::Float64, T::Float64)
    arg = -2 * cos2pi(u_flip - s1) * cos2pi(u_flip - s2) / T
    return 1 - exp(arg)
end

"""
    simulate_xy_wolff(N::Int64, T::Float64, esteps::Inf64, nsteps::Int64; from_infinity=false)

Simulate the XY model for a `(N, N)` lattice at temperature `T` with `eqsteps` number of steps for equilibration and `nsteps` number of steps for measurements.
"""
function simulate_xy_wolff(N::Int64, T::Float64, esteps::Int64, nsteps::Int64; from_infinity = false)
    if from_infinity
        spins = rand(0.0:0.1:1.0, (N, N))
    else
        spins = fill(0.0, (N, N))
    end

    for i = 1:esteps
        xy_wolff_step!(spins, N, T)
    end

    E_arr = zeros(Float64, nsteps)
    for i = 1:nsteps
        xy_wolff_step!(spins, N, T)
        E_arr[i] = total_energy(spins, N)
    end

    return E_arr
end

"""
    plot_spins(spins::Matrix, N::Int64)

Plot `(N, N)` lattice `spins` of XY model as vectors.
"""
function plot_spins(spins::Matrix, N::Int64)
    x = repeat(collect(1:N), outer = N)
    y = repeat(collect(1:N), inner = N)
    vx = zeros(Float64, size(x))
    vy = zeros(Float64, size(x))
    # rotation = zeros(Float64, size(x))
    @inbounds for i in eachindex(x, y)
        vx[i] = cos2pi(spins[x[i], y[i]])
        vy[i] = sin2pi(spins[x[i], y[i]])
        # rotation[i] = mod(abs(0.5 - spins[x[i], y[i]]), 1)
    end

    return arrows(
        x, y, vx, vy,
        arrowsize = 10, lengthscale = 0.5,
        align = :center, normalize = true,
        # arrowcolor=rotation, linecolor=rotation
    )
end

"""
    total_energy(spins)

Calculate the total energy of the square spin lattice `spins` of size `(N, N)` with zero field and J=1.
"""
function total_energy(spins::Matrix, N::Int64)
    running_sum = 0
    for idx in CartesianIndices(spins)
        s_k = spins[idx]
        for δ ∈ CartesianIndex.([(1, 0), (N - 1, 0), (0, 1), (0, N - 1)])
            nn = idx + δ
            nn = CartesianIndex(mod1.(Tuple(nn), N))  # Apply periodic boundary conditions
            running_sum += cos2pi(s_k - spins[nn])
        end
    end
    return -running_sum / 2  # divide by 2 because each bond counted twice
end

"""
    specific_heat(u_vals, T, N)

Calculate the specific heat from given array of internal energy per site (`N²` sites) at temperature `T`.
"""
function specific_heat(E_vals, T, N)
    return (mean(E_vals .^ 2) - mean(E_vals)^2) / (N * T)^2
end

"""
    bootstrap_err(samples, calc_qty; r=100)

Estimate the error in the given samples by bootstrap method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `r` is a keyword arguments giving number of resamples.
"""
function bootstrap_err(samples, calc_qty, args...; r = 100)
    nob = length(samples)
    resample_arr = zeros(Float64, nob)
    for i = 1:r
        resample = rand(samples, nob)
        resample_arr[i] = calc_qty(resample, args...)
    end
    err = std(resample_arr, corrected = false)
    return err
end

"""
    blocking_err(samples, calc_qty; blocks=10)

Estimate the error in the given samples by blocking method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `blocks` is a keyword arguments giving number of blocks.
"""
function blocking_err(samples, calc_qty, args...; blocks = 20)
    block_array = zeros(Float64, blocks)
    blocklength = length(samples) ÷ blocks
    for i = 1:blocks
        sample_block = samples[(i-1)*blocklength+1:i*blocklength]
        block_array[i] = calc_qty(sample_block, args...)
    end
    err = std(block_array)
    return err
end

function cos2pi(x)
    return cos(2 * pi * x)
end

function sin2pi(x)
    return sin(2 * pi * x)
end