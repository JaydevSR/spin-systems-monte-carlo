include("../src/spinmc.jl")

#=
Perform simulation
=#

N = 50
println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

esteps = 2000  # Number of steps for equilibration
nsteps = 10000  # Number of steps for measurements

T = 1.5
println("Calculating spin-spin correlation function for T = $(T) ...")
ss_corrs = zeros(Float64, N)
nsamples = 0

spins = rand(0.0:0.1:1.0, (N, N))

for step = 1:esteps
    xywolff_step!(spins, N, T)
end

for step = 1:nsteps
    global nsamples, ss_corrs
    if step%100 == 0
        println("Step: $(step)/$(nsteps) ...")
    end
    xywolff_step!(spins, N, T)
    if step%100 == 0
        ss_corrs = ss_corrs .+ ss_correlation_fn(spins, N; metric=xy_spindot)
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

save("simulations/results/xy/ss_corr_N$(N)_Temp($(T)).png", f)

println("Program Finished!")
println("===========================\n")