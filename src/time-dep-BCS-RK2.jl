using LinearAlgebra, Plots
include("param.jl")
include("output.jl")

function time_dep_BCS_RK2!(param::Param, output::Output)
    # gap function
    Δ = zeros(ComplexF64, param.Nt+1)
    Δ[1] = 1.0

    # pseudospin
    σk = zeros(Float64, 3, param.Nk, param.Nt+1)

    # Pauli matrices
    τ1 = [0.0 1.0; 1.0  0.0]
    τ2 = [0.0 -im; im   0.0]
    τ3 = [1.0 0.0; 0.0 -1.0]
    τ = [τ1, τ2, τ3]
    
    # BdG Hamiltonian
    function hk(k::Float64, Δ::ComplexF64)
        return (-2.0*param.t_hop*cos(k)-param.μ)*τ3 - Δ*(τ1-im*τ2)/2.0 - conj(Δ)*(τ1+im*τ2)/2.0
    end

    # Effective magnetic field
    function bk(k::Float64, Δ::ComplexF64)
        return [-real(Δ), -imag(Δ), -2.0*param.t_hop*cos(k)-param.μ]
    end

    # equilibrium mean-field self-consistency
    diff = 1.0
    iter = 0
    while diff > param.tolerance
        Δ_temp = 0.0
        for ik in 1:param.Nk
            k = ik*param.dk
            Ek,Uk = eigen(hk(k,Δ[1]))
            f = diagm([1.0/(exp(Ek[1]/param.T)+1.0), 1.0/(exp(Ek[2]/param.T)+1.0)])
            Δ_temp += tr(f*Uk'*(τ1+im*τ2)/2.0*Uk)
        end
        Δ_temp = Δ_temp*param.Vi/param.Nk
        diff = abs(Δ[1]-Δ_temp)
        Δ[1] = Δ_temp
        iter += 1
    end
    output.t_list[1] = 0.0
    output.Δ[1] = Δ[1]
    output.iter_list[1] = iter

    # initial pseudospin configuration
    E_tot = 0.0
    N_tot = 0.0
    for ik in 1:param.Nk
        k = ik*param.dk
        Ek,Uk = eigen(hk(k,Δ[1]))
        f = diagm([1.0/(exp(Ek[1]/param.T)+1.0), 1.0/(exp(Ek[2]/param.T)+1.0)])
        for α in 1:3
            σk[α,ik,1] = real(tr(f*Uk'*τ[α]*Uk)/2.0)
        end
        E_tot += bk(k,Δ[1])⋅σk[:,ik,1]
        N_tot += σk[3,ik,1]
    end
    E_tot = 2.0*E_tot[1]/param.Nk + abs2(Δ[1])/param.Vi
    N_tot = N_tot/param.Nk + 1.0
    output.E_tot[1] = E_tot
    output.N_tot[1] = N_tot

    # time evolution
    for it in 1:param.Nt
        t = it*param.dt
        Δ[it+1] = Δ[it]
        σk[:,:,it+1] = σk[:,:,it]
        diff = 1.0
        iter = 0
        while diff > param.tolerance
            Δ_temp = 0.0
            for ik in 1:param.Nk
                k = ik*param.dk
                # trapezoid on (t, t+Δt)
                σk[:,ik,it+1] = σk[:,ik,it] + param.dt*(bk(k,Δ[it])×σk[:,ik,it]+bk(k,Δ[it+1])×σk[:,ik,it+1])
                Δ_temp += σk[1,ik,it+1] + im*σk[2,ik,it+1]
            end
            Δ_temp = Δ_temp*param.Vf/param.Nk
            diff = abs(Δ[it+1]-Δ_temp)
            Δ[it+1] = Δ_temp
            iter += 1
        end
        output.t_list[it+1] = t
        output.Δ[it+1] = Δ[it+1]
        output.iter_list[it+1] = iter
        E_tot = 0.0
        N_tot = 0.0
        for ik in 1:param.Nk
            k = ik*param.dk
            E_tot += bk(k,Δ[it+1])⋅σk[:,ik,it+1]
            N_tot += σk[3,ik,it+1]
        end
        E_tot = 2.0*E_tot/param.Nk + abs2(Δ[it+1])/param.Vf
        N_tot = N_tot/param.Nk + 1.0
        output.E_tot[it+1] = E_tot
        output.N_tot[it+1] = N_tot
    end
end
