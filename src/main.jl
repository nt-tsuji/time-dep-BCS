include("param.jl")
include("output.jl")
include("time-dep-BCS-RK2.jl") # second-order implicit Runge-Kutta
# include("time-dep-BCS-RK4.jl") # fourth-order implicit Runge-Kutta

Nt = 1000
tmax = 50.0
Vi = 2.0
Vf = 6.0
t_hop = 1.0
T = 0.05
Nk = 200
tolerance = 0.000001

param = Param(Nt, tmax, Vi, Vf, t_hop, T, Nk, tolerance)
output = Output(param.Nt)
time_dep_BCS_RK2!(param, output)
# time_dep_BCS_RK4!(param, output)

fig_Delta = plot(output.t_list, [real(output.Δ), imag(output.Δ)], xlim=(0,param.tmax), ylim=(0,4), xlabel="t", label=["Re Δ" "Im Δ"])
fig_E_tot = plot(output.t_list, output.E_tot, xlim=(0,param.tmax), ylim=(-1.5,-1.3), xlabel="t", label="E_tot")

savefig(fig_Delta,"Delta.png")
savefig(fig_E_tot,"E_tot.png")
