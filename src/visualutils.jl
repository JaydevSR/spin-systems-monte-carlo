"""
    ising_plot_spins(spins::Matrix, N::Int64)

Plot `(N, N)` lattice `spins` of Ising model as binary heatmap.
"""
function xy_plot_spins(spins::Matrix, N::Int64; colorscheme=:winter)
    return heatmap(1:N, 1:N, spins, colormap=cgrad(colorscheme, 2, categorical=true))
end

"""
    xy_plot_spins(spins::Matrix, N::Int64)

Plot `(N, N)` lattice `spins` of XY model as vectors.
"""
function xy_plot_spins(spins::Matrix, N::Int64)
    x = repeat(collect(1:N), outer = N)
    y = repeat(collect(1:N), inner = N)
    vx = zeros(Float64, size(x))
    vy = zeros(Float64, size(x))
    rotation = zeros(Float64, size(x))
    @inbounds for i in eachindex(x, y)
        vx[i] = cos2pi(spins[x[i], y[i]])
        vy[i] = sin2pi(spins[x[i], y[i]])
        rotation[i] = mod(abs(0.5 - spins[x[i], y[i]]), 1)
    end

    return arrows(
        x, y, vx, vy,
        arrowsize = 10, lengthscale = 0.5,
        align = :center, normalize = true,
        arrowcolor=rotation, linecolor=rotation
    )
end