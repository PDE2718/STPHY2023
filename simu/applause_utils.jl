# logistic functions
logistic(x, τ) = 1 / (1 + exp(-x / τ))
logistic(x, τ, x0) = 1 / (1 + exp(-(x-x0) / τ))

# generating gaussian distributed variables

function randnσc(σc)
    r = randn()
    while abs(r) > σc
        r = randn()
    end
    return r
end
function gaussian_param(pc,σ_rel,cut_rel)
    return pc * (1 + σ_rel * randnσc(cut_rel))
end
gaussian_param(pc, σ_rel) = gaussian_param(pc, σ_rel, 1.0)

# window average of signal
function windowavg(t, y, τ)
    @assert length(t) == length(y)
    y1 = deepcopy(y)
    t1 = deepcopy(t)
    y2 = reverse(y)
    t2 = -reverse(t)
    @inbounds for i ∈ 2:lastindex(t1)
        y1[i] = y1[i-1] + (y1[i] - y1[i-1]) * (t1[i] - t1[i-1]) / τ
    end
    for i ∈ 2:lastindex(t2)
        y2[i] = y2[i-1] + (y2[i] - y2[i-1]) * (t2[i] - t2[i-1]) / τ
    end
    return (y1 .+ reverse(y2)) ./ 2
end

function F_indv(u)
    return 2.0logistic(real(u) - 1.0, 0.1)
end

function F_mean(u)
    return mean(F_indv, u)
end

function sync_term(ϕ, u)
    return abs(ϕ) * sin(angle(ϕ / u))
end

function gγ(u, g0, τ, F, t)
    g = g0 * F
    γ = 1 / τ
    p = logistic(t - τ, 1.0)
    bound = logistic(abs(u) - 1, 0.2)
    gγ = (1 - p) * (1 - bound) * g - p * γ
    return gγ
end

function order_param(u)
    return mean(u)/√mean(abs2,u)
end


function applause(du, u, p, t)
    ω, τ, g, κ = p
    F = F_mean(u)
    ϕ = order_param(u)
    for j ∈ eachindex(du, u, ω, τ, g, κ)
        du[j] = 1im * (ω[j] + κ[j] * sync_term(ϕ, u[j])) + gγ(u[j], g[j], τ[j], F, t)
        du[j] *= u[j]
    end
end

