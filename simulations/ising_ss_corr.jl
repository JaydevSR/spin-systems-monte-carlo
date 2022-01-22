include("../src/spinmc.jl")

N = 100  # Lattice size

println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

# Initialize lattice
spins = ones(N, N)  # T = 0
# spins = rand([-1, 1], (N, N))  # T = ∞

eqsteps = 2000  # Number of steps for equilibration
nsteps = 10000  # Number of steps for measurements

T = round(ising_Tc, digits=4)
println("Calculating spin-spin correlation function for T = $(T) ...")
ss_corrs = zeros(Float64, N)
nsamples = 0

E = ising_total_energy(spins)
M = ising_total_magnetization(spins)

# Let the system reach equilibrium
for step = 1:eqsteps
    global E, M
    E, M = isingmetro_step!(spins, T, E, M)
end

# Iterate for calculating averages
for step = 1:nsteps
    global E, M, nsamples, ss_corrs
    if step%100 == 0
        println("Step: $(step)/$(nsteps) ...")
    end
    E, M = isingmetro_step!(spins, T, E, M)
    if step%100 == 0
        ss_corrs = ss_corrs .+ ss_correlation_fn(spins, N)
        nsamples += 1
    end
end

ss_corrs = ss_corrs ./ nsamples
ss_corrs = (ss_corrs[1:end] .+ ss_corrs[end:-1:1]) ./ 2


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "r = |i-j|", ylabel = "<s_i * s_j>",
    title = "Spin-spin correlation function (Lattice Size $(N), T=$(T)")

lines!(ax1, 0:N-1, ss_corrs[1:N])
scatter!(
    ax1, 0:N-1, ss_corrs[1:N],
    markersize = 10
)

save("simulations/results/ising/ss_corr_N$(N)_Temp($(T)).png", f)

println("Program Finished!")
println("===========================\n")