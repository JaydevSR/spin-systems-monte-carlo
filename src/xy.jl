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
    while !isempty(stack)
        k = pop!(stack)
        kval = spins[k...]
        xywolff_flip_spin!(spins, k, u_flip)
        @inbounds for δ ∈ ([1, 0], [N - 1, 0], [0, 1], [0, N - 1])
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            nnval = spins[nn...]
            if !cluster[nn...] && rand() < xywolff_Padd(u_flip, nnval, kval, T)
                push!(stack, nn)
                cluster[nn...] = true
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
    new = mod(new + 1, 1)
    spins[pos...] = new
    return old, spins[pos...]
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
    x = along_u .* cos2pi(u_flip) .- normal_u .* sin2pi(u_flip)
    y = along_u .* sin2pi(u_flip) .+ normal_u .* cos2pi(u_flip)
    XYlattice = mod1.(atan.(y, x) ./ (2pi), 1)
    return XYlattice
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