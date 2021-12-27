using CairoMakie

function xy_wolff_step!(spins, N, T)
    seed = rand(1:N, 2)  # seed spin position
    u_flip = rand()  # Random unit vector in xy plane
    cluster = get_cluster!(spins, seed, u_flip, T)
end

function get_cluster!(spins, seed, u_flip, T)
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
    return cluster
end

function P_add(u_flip, s1, s2, T)
    arg = -2 * cos(u_flip - s1) * cos(u_flip - s2) / T
    return 1 - exp(arg)
end

function flip_spin!(spins, pos, u_flip)
    old = spins[pos...]
    new = 0.5 + 2*u_flip - old  # flipping w.r.t vector with angle ϕ: θ --> π + 2ϕ - θ
    new = mod(new, 1)
    spins[pos...] = new
    return old, spins[pos...]
end

function simulate_xy_wolff(N, T, esteps, nsteps; from_infinity=false)
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
        # u_arr[i] = total_energy(spins) / N^2
    end

    return u_arr
end

function plot_spins(spins, N)
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

# simulation

N = 20
T = 0.8

# u_arr = simulate_xy_wolff(N, T, 2000, 6000; from_infinity=true)
spins = rand(Float64, (N, N))

# h = heatmap(spins, c=:twilight)
h = plot_spins(spins, N)
display(h)

for i=1:6000
    if i%1000==0
        h = plot_spins(spins, N)
        display(h)
    end
    xy_wolff_step!(spins, N, T)
end