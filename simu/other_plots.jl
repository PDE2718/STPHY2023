using Plots, LaTeXStrings, Measures
default(framestyle = :box, labels = nothing)

begin
    plt4 = plot(
        xlabel = L"t"* " (s)",
        ylims=(0.0, 2.6),
        ylabel="individual clapping amplitude",
        xlims=(0,10),
        size = (700,300),
        margins = 3mm,
    )
    plot!(u_t, [F_indv.(sol[i, :]) for i ∈ 1:4], fill = 0, labels = "person ".*["$(i)" for i ∈ 1:4] |> permutedims)
end
savefig(plt4, "individualclap.pdf")

ttg = LinRange(0, 10, 101)
FFg = LinRange(0, 0.4, 100)
g_mesh = [gγ(complex(1.1), 4.0, 5.0, F, t) for F ∈ FFg, t ∈ ttg]
plt5 = surface(ttg, FFg, g_mesh,
    xlabel = L"t"* " (s)",
    ylabel= L"\overline{F}",
    zlabel= L"g_j(\overline{F},u_j,t)",
    size=(600, 400),
)
savefig(plt5, "gfunc.pdf")

begin
    ωdist = [gaussian_param(3.0, 0.2) for i ∈ 1:400]
end

plt_ωdist = plot(
    xlabel=L"\omega_j/2\pi" * " (Hz)",
    ylims=(0.0, 1.0),
    ylabel="PDF",
    xlims=(2, 4),
    size=(350, 250),
    margins=0mm,
)
histogram!(ωdist, normalize=true, bins=:fd)

savefig(plt_ωdist,"omegadist.pdf")