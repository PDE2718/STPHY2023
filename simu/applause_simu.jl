cd(@__DIR__);
using LinearAlgebra                             # math
using DifferentialEquations.OrdinaryDiffEq      # ODEs
using Statistics, StatsBase, FFTW               # Statistics & FFT
using Plots, LaTeXStrings, Measures             # Plots
default(framestyle=:box, labels=nothing)

include("applause_utils.jl")

function simu_applause(para; tspan=(0.0, 15.0), dt=0.02)
    N = para[:N]
    σ0 = para[:σ0]
    ω = [gaussian_param(para[:ω0], σ0) for _ ∈ 1:N]
    τ = [gaussian_param(para[:τ0], σ0) for _ ∈ 1:N]
    g = [gaussian_param(para[:g0], σ0) for _ ∈ 1:N]
    κ = [gaussian_param(para[:κ0], σ0) for _ ∈ 1:N]
    u0 = [0.5cis(2π * rand()) for i ∈ 1:N]
    for i ∈ eachindex(u0)
        u0[i] *= (rand() < para[:η0] ? 2.0 : 1.0)
    end
    p = (ω, τ, g, κ)
    prob = ODEProblem(applause, u0, tspan, p)
    sol = solve(prob, Tsit5(); dt=dt, adaptive=false, dense=false)
    return sol
end

begin
    para = (N=400, σ0=0.3, τ0=5.0, g0=4.0, κ0=5.30, η0=0.5, ω0=2π * 3)
    sol = simu_applause(para)
    tspan = (u_t[1], u_t[end])
    u_t = sol.t
    u_F = F_mean.(sol.u)
    u_ϕ = order_param.(sol.u)
    u_Fs = windowavg(u_t, u_F, 0.8)
end

begin
plt1 = plot(
    ylims=(0.0, 0.8),
    ylabel="sound intensity  " * L"\bar{F}(t)",
    xlims=tspan,
    label=L"F(t)",
)
plot!(u_t, u_Fs,
    lw=2, lc=4, fill=0, fillcolor=4, fillalpha=0.2,
    label="filtered  " * L"T_w=0.8",
)
plot!(u_t, u_F,
    lw=1.5, lc=1,
    label=L"\bar{F}(t)",
)
end

function stft(yg, tg, t, σ)
    yp = @. yg * exp(-0.5 * ((tg - t) / σ)^2) / (σ * √(2π))
    return abs2.(fft(yp))
end
fs = 1 / (u_t[2] - u_t[1])
tg_tf = LinRange(tspan..., 200)
fg = LinRange(0, fs, length(u_t))
tf_spec = hcat([stft(u_F, u_t, t, 1.0) for t ∈ tg_tf]...)

plt2 = heatmap(tg_tf, fg, log.(tf_spec),
    ylims=(0, 6.0), clims=(-8, 2), xlims=(0, 13.2),
    xguide=L"t" * " (s)", yguide=L"f" * " (Hz)",
    colorbar_title="spectrum intensity " * L"\log[I(f)]",
)

plt0 = plot(plt1, plt2, layout=grid(2, 1), size=(700, 500), show=true, margins=1mm, right_margin=2mm)

# savefig(plt0, "xxx.pdf")