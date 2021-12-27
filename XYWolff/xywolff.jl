using CairoMakie
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
        for δ ∈ ([1, 0], [N-1, 0], [0, 1], [0, N-1])
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            nnval = spins[nn...]
            if abs(sval - nnval) <= 0.25 && !cluster[nn...] && rand() < P_add(u_flip, nnval, sval, T)
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
    arg = -2 * cos(u_flip - s1) * cos(u_flip - s2) / T
    return 1 - exp(arg)
end

"""
    flip_spin!(spins::Matrix, pos::AbstractArray, u_flip::Float64)

Flip the spin at position `pos` inside lattice `spins` w.r.t. angle `u_flip`.
"""
function flip_spin!(spins::Matrix, pos::AbstractArray, u_flip::Float64)
    old = spins[pos...]
    new = 0.5 + 2*u_flip - old  # flipping w.r.t vector with angle ϕ: θ --> π + 2ϕ - θ
    new = mod(new, 1)
    spins[pos...] = new
    return old, spins[pos...]
end

"""
    simulate_xy_wolff(N::Int64, T::Float64, esteps::Inf64, nsteps::Int64; from_infinity=false)

Simulate the XY model for a `(N, N)` lattice at temperature `T` with `eqsteps` number of steps for equilibration and `nsteps` number of steps for measurements.
"""
function simulate_xy_wolff(N::Int64, T::Float64, esteps::Int64, nsteps::Int64; from_infinity=false)
    if from_infinity
        spins = rand(Float64, (N, N))
    else
        spins = zeros(Float64, (N, N))
    end

    for i=1:esteps
        xy_wolff_step!(spins, N, T)
    end

    u_arr = zeros(Float64, nsteps)
    for i=1:nsteps
        xy_wolff_step!(spins, N, T)
        u_arr[i] = total_energy(spins, N) / N^2
    end

    return u_arr
end

"""
    plot_spins(spins::Matrix, N::Int64)

Plot `(N, N)` lattice `spins` of XY model as vectors.
"""
function plot_spins(spins::Matrix, N::Int64)
    x = repeat(collect(1:N), outer=N)
    y = repeat(collect(1:N), inner=N)
    vx = zeros(Float64, size(x))
    vy = zeros(Float64, size(x))
    # rotation = zeros(Float64, size(x))
    @inbounds for i in eachindex(x, y)
        vx[i] = cos(2*pi*spins[x[i], y[i]])
        vy[i] = sin(2*pi*spins[x[i], y[i]])
        # rotation[i] = mod(abs(0.5 - spins[x[i], y[i]]), 1)
    end
    
    return arrows(
        x, y, vx, vy, 
        arrowsize = 10, lengthscale = 0.5, 
        align=:center, normalize=true,
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
        for δ ∈ CartesianIndex.([(1, 0), (N-1, 0), (0, 1), (0, N-1)])
            nn = idx + δ
            nn = CartesianIndex(mod1.(Tuple(nn), N))  # Apply periodic boundary conditions
            running_sum += cos(s_k-spins[nn])
        end
    end
    return -running_sum / 2  # divide by 2 because each bond counted twice
end

"""
    specific_heat(u_vals, T, N)

Calculate the specific heat from given array of internal energy per site (`N²` sites) at temperature `T`.
"""
function specific_heat(u_vals, T, N)
    return (T^-2) * N^2 * var(u_vals, corrected=false)
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

# simulation

N = 10

Temps = [i for i=0.1:0.1:2.4]
esteps = 2000  # Number of steps for equilibration
nsteps = 4000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

c_T = zeros(Float64, length(Temps))  # Array of specific heat
err_c_T = zeros(Float64, length(Temps))

for i=1:length(Temps)
    global spins
    T = Temps[i]
    
    u_arr = simulate_xy_wolff(N, T, esteps, nsteps; from_infinity=true)

    u_T[i] = mean(u_arr)
    err_u_T[i] = blocking_err(u_arr, mean)

    c_T[i] = specific_heat(u_arr, T, N)
    err_c_T[i] = blocking_err(u_arr, specific_heat, T, N)
end


#=
Plots
=#
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "XY Model for Lattice Size $(N)")

ax2 = Axis(f[1, 2], xlabel = "temperature, T", ylabel = "specific heat, c",
    title = "XY Model for Lattice Size $(N)")

errorbars!(ax1, Temps, u_T, err_u_T)
scatter!(ax1, Temps, u_T, yerr=err_u_T)
errorbars!(ax2, Temps, c_T, err_c_T)
scatter!(ax2, Temps, c_T, yerr=err_c_T)

save("XYWolff/plots/u&c_vs_T_$(N).png", f)