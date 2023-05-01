using AbstractGPs
using StatsPlots, StatsPlots.PlotMeasures
using Random

Random.seed!(2)

f = GP(Matern32Kernel())
xx = range(0, 10, length=100)
fx = f(xx, 0.01)

yl = (-3, 3)
size = (800, 300)
p = plot(xx, rand(fx), legend=false, margin=15px, size=size, ylims=yl, xlabel="X", ylabel="Y")
savefig(joinpath(@__DIR__, "unconditional_1.png"))
for i in 1:5
    plot!(p, xx, rand(fx), color=1, ylims=yl)
end
p
savefig(joinpath(@__DIR__, "unconditional.png"))


x = range(0, 10, length=5)
y = [-0.5, 0.1, -1.4, 1, 1.1]

scatter(x, y, ylims=yl, xlabel="X", ylabel="Y", color=2, label="",
    size=size, margin=15px)
savefig(joinpath(@__DIR__, "data.png"))

p_fx = posterior(f(x, 1e-2), y)
plot(xx, p_fx, label="",
    size=size, margin=15px)
scatter!(x, y, ylims=yl, xlabel="X", ylabel="Y", color=2, label="")
savefig(joinpath(@__DIR__, "kriging.png"))

p2 = plot(xx, rand(p_fx(xx, 0.01)), color=1, alpha=0.5, legend=false,
    ylims=yl, xlabel="X", ylabel="Y", size=size, margin=15px)
for i in 1:15
    plot!(p2, xx, rand(p_fx(xx, 0.01)), color=1, alpha=0.5)
end
scatter!(p2, x, y, color=2)
savefig(joinpath(@__DIR__, "conditional.png"))


