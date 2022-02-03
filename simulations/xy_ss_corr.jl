include("../src/spinmc.jl")

#=
Perform simulation
=#

N = 50
println("================================\n")
println("    Lattice Size: $(N) x $(N)")
println("================================\n")

# T > Tc
# T = 1.5
# esteps = 4000
# nsteps = 40000

# T < Tc
T = 0.4
esteps = 1000  # Number of steps for equilibration
nsteps = 10000  # Number of steps for measurements

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
# @. ss_corrs[2:end] = (ss_corrs[2:end] + ss_corrs[end:-1:2]) / 2

#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "r = |i-j|", ylabel = "<s_i * s_j>",
    yminorticksvisible = true, yminorgridvisible = true,
    title = "Spin-spin correlation function (Lattice Size $(N), T=$(T))"
)

ax2 = Axis(f[2, 1], ylabel = "ln(<s_i * s_j>)", yscale=log, 
    # xlabel = "r = |i-j|",  # For exponential decay
    xlabel = "ln(r) = ln(|i-j|)", xscale=log,  # For power law decay
    yminorticksvisible = true, yminorgridvisible = true,
    # title = "xscale=linear, yscale=logscale"  # For exponential decay
    title = "xscale=logscale, yscale=logscale"  # For power law decay
)

scatter!(
        ax1, 0:N÷2-1, ss_corrs[1:N÷2],
        markersize = 7,
)

lines!(ax2, 1:N÷2-1, ss_corrs[2:N÷2])
scatter!(
    ax2, 1:N÷2-1, ss_corrs[2:N÷2],
    markersize = 7
)

save("simulations/results/xy/ss_corr_N$(N)_Temp($(T)).png", f)

println("Program Finished!")
println("===========================\n")