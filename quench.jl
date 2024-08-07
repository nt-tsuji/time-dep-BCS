using LinearAlgebra, Plots

# parameters
struct Param
    Nt::Int64
    tmax::Float64
    Vi::Float64
    Vf::Float64
    t_hop::Float64
    T::Float64
    Nk::Int64
    tolerance::Float64
end

mutable struct Output
    t_list::Array{Float64,1}
    Δ_re::Array{Float64,1}
    Δ_im::Array{Float64,1}
    E_tot::Array{Float64,1}
    iter_list::Array{Int16,1}
    function Output(Nt::Int64)
        new(zeros(Float64,Nt+1), zeros(Float64,Nt+1), zeros(Float64,Nt+1), zeros(Float64,Nt+1), zeros(Int64,Nt+1))
    end
end

function time_dep_BCS(param::Param, output::Output)
    dk = 2*π/param.Nk
    dt = param.tmax/param.Nt

    # gap function
    Δ = zeros(ComplexF64, param.Nt+1)
    Δ[1] = 1.

    # pseudospin
    σk = zeros(Float64, 3, param.Nk, param.Nt+1)

    # Pauli matrices
    τ1 = [0.  1.; 1.  0.]
    τ2 = [0. -im; im  0.]
    τ3 = [1.  0.; 0. -1.]
    τ = [τ1, τ2, τ3]
    
    # BdG Hamiltonian
    function hk(k::Float64, Δ::ComplexF64)
        return -2*param.t_hop*cos(k)*τ3 - Δ*(τ1-im*τ2)/2 - conj(Δ)*(τ1+im*τ2)/2
    end

    # Effective magnetic field
    function bk(k::Float64, Δ::ComplexF64)
        return [-real(Δ), -imag(Δ), -2*param.t_hop*cos(k)]
    end

    # equilibrium mean-field self-consistency
    diff = 1.
    iter = 0
    while diff > param.tolerance
        Δ_temp = 0.
        for ik in 1:param.Nk
            k = ik*dk
            Ek,Uk = eigen(hk(k,Δ[1]))
            f = diagm([1/(exp(Ek[1]/param.T)+1.0), 1/(exp(Ek[2]/param.T)+1.0)])
            Δ_temp += tr(f*Uk'*(τ1+im*τ2)/2*Uk)
        end
        Δ_temp = Δ_temp*param.Vi/param.Nk
        diff = abs(Δ[1]-Δ_temp)
        Δ[1] = Δ_temp
        iter += 1
    end
    output.t_list[1] = 0.
    output.Δ_re[1] = real(Δ[1])
    output.Δ_im[1] = imag(Δ[1])
    output.iter_list[1] = iter

    # initial pseudospin configuration
    E_tot = 0.
    for ik in 1:param.Nk
        k = ik*dk;
        Ek,Uk = eigen(hk(k,Δ[1]))
        f = diagm([1/(exp(Ek[1]/param.T)+1.0), 1/(exp(Ek[2]/param.T)+1.0)])
        for α in 1:3
            σk[α,ik,1] = real(tr(f*Uk'*τ[α]*Uk)/2)
        end
        E_tot += bk(k,Δ[1])⋅σk[:,ik,1]
    end
    E_tot = 2*E_tot[1]/param.Nk + abs2(Δ[1])/param.Vi
    output.E_tot[1] = E_tot

    # time evolution
    for it in 1:param.Nt
        t = it*dt
        Δ[it+1] = Δ[it]
        σk[:,:,it+1] = σk[:,:,it]
        diff = 1.
        iter = 0
        while diff > param.tolerance
            Δ_temp = 0.
            for ik in 1:param.Nk
                k = ik*dk
                σk[:,ik,it+1] = σk[:,ik,it] + dt*(bk(k,Δ[it])×σk[:,ik,it]+bk(k,Δ[it+1])×σk[:,ik,it+1])
                Δ_temp += σk[1,ik,it+1] + im*σk[2,ik,it+1]
            end
            Δ_temp = Δ_temp*param.Vf/param.Nk
            diff = abs(Δ[it+1]-Δ_temp)
            Δ[it+1] = Δ_temp
            iter += 1
        end
        output.t_list[it+1] = t
        output.Δ_re[it+1] = real(Δ[it+1])
        output.Δ_im[it+1] = imag(Δ[it+1])
        output.iter_list[it+1] = iter
        E_tot = 0.
        for ik in 1:param.Nk
            k = ik*dk;
            E_tot += bk(k,Δ[it+1])⋅σk[:,ik,it+1]
        end
        E_tot = 2*E_tot/param.Nk + abs2(Δ[it+1])/param.Vf
        output.E_tot[it+1] = E_tot
    end
end

Nt = 1000;
tmax = 50.;
Vi = 2.;
Vf = 6.;
t_hop = 1.;
T = 0.05;
Nk = 200;
tolerance = 0.000001;

param = Param(Nt, tmax, Vi, Vf, t_hop, T, Nk, tolerance);
output = Output(param.Nt);
time_dep_BCS(param, output);

fig_Delta = plot(output.t_list, [output.Δ_re, output.Δ_im], xlim=(0,param.tmax), ylim=(0,4), xlabel="t", label=["Re Δ" "Im Δ"])
fig_E_tot = plot(output.t_list, output.E_tot, xlim=(0,param.tmax), ylim=(-1.5,-1.3), xlabel="t", label="E_tot")

savefig(fig_Delta,"Delta.png")
savefig(fig_E_tot,"E_tot.png")
