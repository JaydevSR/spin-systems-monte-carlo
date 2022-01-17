using Statistics


"""
    xy_wolff_step!(spins::Matrix, N::Int64, T::Float64)

Perform one step of Wolff algorithm for XY model (lattice `spins` of size `(N, N)` at temperature `T`).
"""
function xy_wolff_step!(spins::Matrix, N::Int64, T::Float64)
    seed = rand(1:N, 2)  # seed spin position
    u_flip = rand()  # Random unit vector in xy plane
    cluster_update!(spins, seed, u_flip, T)
end


"""
    cluster_update!(spins::Matrix, seed::AbstractArray, u_flip::Float64, T::Float64)

Build a cluster among `spins` starting at `seed` at temperature `T`. Flip the cluster w.r.t angle `u_flip`.
"""
function cluster_update!(spins::Matrix, seed::AbstractArray, u_flip::Float64, T::Float64)
    cluster = falses(size(spins))
    sval = spins[seed...]
    stack = [seed]
    cluster[seed...] = true
    flip_spin!(spins, seed, u_flip)
    while !isempty(stack)
        k = pop!(stack)
        @inbounds for δ ∈ ([1, 0], [N - 1, 0], [0, 1], [0, N - 1])
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            nnval = spins[nn...]
            if isparallel(nnval, sval) && !cluster[nn...] && rand() < P_add(u_flip, nnval, sval, T)
                push!(stack, nn)
                cluster[nn...] = true
                flip_spin!(spins, nn, u_flip)
            end
        end
    end
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
    flip_spin!(spins::Matrix, pos::AbstractArray, u_flip::Float64)

Flip the spin at position `pos` inside lattice `spins` w.r.t. angle `u_flip`.
"""
function flip_spin!(spins::Matrix, pos::AbstractArray, u_flip::Float64)
    old = spins[pos...]
    new = 0.5 + 2 * u_flip - old  # flipping w.r.t vector with angle ϕ: θ --> π + 2ϕ - θ
    new = mod1(new, 1)
    spins[pos...] = new
    return old, spins[pos...]
end

function isparallel(s1, s2)
    tht1, tht2 = 2pi*s1, 2pi*s2
    if cos(tht1 - tht2) > 0
        return true
    end
    return false
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