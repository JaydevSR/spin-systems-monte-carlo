using CairoMakie
using Statistics
using SpecialFunctions

include("ising.jl")
include("xy.jl")
include("observables.jl")
include("statutils.jl")
include("visualutils.jl")
include("isingexact.jl")

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
        xywolff_step!(spins, N, T)
    end

    E_arr = zeros(Float64, nsteps)
    for i = 1:nsteps
        xywolff_step!(spins, N, T)
        E_arr[i] = xy_total_energy(spins, N)
    end

    return E_arr
end