"""
    xywolff_step!(spins::Matrix, N::Int64, T::Float64)

Perform one step of Wolff algorithm for XY model (lattice `spins` of size `(N, N)` at temperature `T`).
"""
function xywolff_step!(spins::Matrix, N::Int64, T::Float64)
    seed = rand(1:N, 2)  # seed spin position
    u_flip = rand()  # Random unit vector in xy plane
    xywolff_cluster_update!(spins, seed, u_flip, T)
end

"""
    xywolff_cluster_update!(spins::Matrix, seed::AbstractArray, u_flip::Float64, T::Float64)

Build a cluster among `spins` starting at `seed` at temperature `T`. Flip the cluster w.r.t angle `u_flip`.
"""
function xywolff_cluster_update!(spins::Matrix, seed::AbstractArray, u_flip::Float64, T::Float64)
    cluster = falses(size(spins))
    sval = spins[seed...]
    stack = [seed]
    cluster[seed...] = true
    xywolff_flip_spin!(spins, seed, u_flip)
    while !isempty(stack)
        k = pop!(stack)
        kval = spins[k...]
        @inbounds for δ ∈ ([1, 0], [N - 1, 0], [0, 1], [0, N - 1])
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            nnval = spins[nn...]
            if !cluster[nn...] && rand() < xywolff_Padd(u_flip, nnval, kval, T)
                push!(stack, nn)
                cluster[nn...] = true
                xywolff_flip_spin!(spins, nn, u_flip)
            end
        end
    end
end

"""
    xywolff_flip_spin!(spins::Matrix, pos::AbstractArray, u_flip::Float64)

Flip the spin at position `pos` inside lattice `spins` w.r.t. angle `u_flip`.
"""
function xywolff_flip_spin!(spins::Matrix, pos::AbstractArray, u_flip::Float64)
    old = spins[pos...]
    new = 0.5 + 2 * u_flip - old  # flipping w.r.t vector with angle ϕ: θ --> π + 2ϕ - θ
    new = mod1(new, 1)
    spins[pos...] = new
    return old, spins[pos...]
end

"""
    xywolff_Padd(u_flip::Float64, s1::Float64, s2::Float64, T::Float64; J=1)

Calculate the probability of adding spin `s2` to cluster of as a neighbour of `s1` at temperature `T` w.r.t angle `u_flip`. The interaction energy defaults to `1`.
"""
function xywolff_Padd(u_flip::Float64, s1::Float64, s2::Float64, T::Float64; J=1)
    arg = -2 * J * cos2pi(u_flip - s1) * cos2pi(u_flip - s2) / T
    return 1 - exp(arg)
end

"""
    xywolff_isparallel(s1, s2)

Determine whether `s1` and `s2` have a component along the same direction or not. Returns `true` or `false`.
"""
function xywolff_isparallel(s1::Float64, s2::Float64)
    tht1, tht2 = 2pi*s1, 2pi*s2
    if 0 <= mod1(tht1 - tht2, 2pi) < pi
        return true
    end
    return false
end

function cos2pi(x::Float64)
    return cos(2 * pi * x)
end

function sin2pi(x::Float64)
    return sin(2 * pi * x)
end

function xy_spindot(s1, s2)
    return cos2pi(s1-s2)
end