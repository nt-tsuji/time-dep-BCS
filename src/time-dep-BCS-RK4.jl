using LinearAlgebra, Plots
include("param.jl")
include("output.jl")

function time_dep_BCS_RK4!(param::Param, output::Output)
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

    # time evolution from t=0 to 2*Δt
    Δ[2] = Δ[1]
    Δ[3] = Δ[1]
    σk[:,:,2] = σk[:,:,1]
    σk[:,:,3] = σk[:,:,1]
    diff = 1.0
    iter = 0
    while diff > param.tolerance
        Δ_temp2 = 0.0
        Δ_temp3 = 0.0
        # Δ(Δt/2)
        Δ_half = 3.0/8.0*Δ[1] + 3.0/4.0*Δ[2] - 1.0/8.0*Δ[3]
        for ik in 1:param.Nk
            k = ik*param.dk
            # σk(Δt/2)
            σk_half = zeros(Float64,3)
            σk_half[:] = 3.0/8.0*σk[:,ik,1] + 3.0/4.0*σk[:,ik,2] - 1.0/8.0*σk[:,ik,3]
            # Simpson on (0, Δt)
            σk[:,ik,2] = σk[:,ik,1] + 2.0*param.dt/6.0*(bk(k,Δ[1])×σk[:,ik,1]+4.0*bk(k,Δ_half)×σk_half[:]+bk(k,Δ[2])×σk[:,ik,2])
            # Simpson on (0, 2*Δt)
            σk[:,ik,3] = σk[:,ik,1] + 2.0*param.dt/3.0*(bk(k,Δ[1])×σk[:,ik,1]+4.0*bk(k,Δ[2])×σk[:,ik,2]+bk(k,Δ[3])×σk[:,ik,3])
            Δ_temp2 += σk[1,ik,2] + im*σk[2,ik,2]
            Δ_temp3 += σk[1,ik,3] + im*σk[2,ik,3]
        end
        Δ_temp2 = Δ_temp2*param.Vf/param.Nk
        Δ_temp3 = Δ_temp3*param.Vf/param.Nk
        diff = abs(Δ[2]-Δ_temp2) + abs(Δ[3]-Δ_temp3)
        Δ[2] = Δ_temp2
        Δ[3] = Δ_temp3
        iter += 1
    end
    output.t_list[2] = param.dt
    output.t_list[3] = 2.0*param.dt
    output.Δ[2] = Δ[2]
    output.Δ[3] = Δ[3]
    output.iter_list[2] = iter
    output.iter_list[3] = iter
    E_tot2 = 0.0
    E_tot3 = 0.0
    N_tot2 = 0.0
    N_tot3 = 0.0
    for ik in 1:param.Nk
        k = ik*param.dk
        E_tot2 += bk(k,Δ[2])⋅σk[:,ik,2]
        E_tot3 += bk(k,Δ[3])⋅σk[:,ik,3]
        N_tot2 += σk[3,ik,2]
        N_tot3 += σk[3,ik,3]
    end
    E_tot2 = 2.0*E_tot2/param.Nk + abs2(Δ[2])/param.Vf
    E_tot3 = 2.0*E_tot3/param.Nk + abs2(Δ[3])/param.Vf
    N_tot2 = N_tot2/param.Nk + 1.0
    N_tot3 = N_tot3/param.Nk + 1.0
    output.E_tot[2] = E_tot2
    output.E_tot[3] = E_tot3
    output.N_tot[2] = N_tot2
    output.N_tot[3] = N_tot3
    
    # time evolution from t to t+Δt
    for it in 3:param.Nt
        t = it*param.dt
        Δ[it+1] = Δ[it]
        σk[:,:,it+1] = σk[:,:,it]
        diff = 1.0
        iter = 0
        while diff > param.tolerance
            Δ_temp = 0.0
            for ik in 1:param.Nk
                k = ik*param.dk
                # Simpson 3/8 on (0, 3*Δt) - Simpson on (0, 2*Δt)
                if (it==3)
                    σk[:,ik,it+1] = σk[:,ik,it] + 2.0*param.dt*(
                        1.0/24.0*bk(k,Δ[1])×σk[:,ik,1]
                        -5.0/24.0*bk(k,Δ[2])×σk[:,ik,2]
                        +19.0/24.0*bk(k,Δ[3])×σk[:,ik,3]
                        +3.0/8.0*bk(k,Δ[4])×σk[:,ik,4])
                # composite Simpson on (0, 4*Δt) - Simpson 3/8 on (0, 3*Δt)
                elseif (it==4)
                    σk[:,ik,it+1] = σk[:,ik,it] + 2.0*param.dt*(
                        -1.0/24.0*bk(k,Δ[1])×σk[:,ik,1]
                        +5.0/24.0*bk(k,Δ[2])×σk[:,ik,2]
                        -11.0/24.0*bk(k,Δ[3])×σk[:,ik,3]
                        +23.0/24.0*bk(k,Δ[4])×σk[:,ik,4]
                        +1.0/3.0*bk(k,Δ[5])×σk[:,ik,5])
                # Gregory on (0, 5*Δt) - composite Simpson on (0, 4*Δt)
                elseif (it==5)
                    σk[:,ik,it+1] = σk[:,ik,it] + 2.0*param.dt*(
                        1.0/24.0*bk(k,Δ[1])×σk[:,ik,1]
                        -1.0/6.0*bk(k,Δ[2])×σk[:,ik,2]
                        +7.0/24.0*bk(k,Δ[3])×σk[:,ik,3]
                        -3.0/8.0*bk(k,Δ[4])×σk[:,ik,4]
                        +5.0/6.0*bk(k,Δ[5])×σk[:,ik,5]
                        +3.0/8.0*bk(k,Δ[6])×σk[:,ik,6])
                # Gregory on (0, t+Δt) - Gregory on (0, t)
                elseif (it>=6)
                    σk[:,ik,it+1] = σk[:,ik,it] + 2.0*param.dt*(
                        1.0/24.0*bk(k,Δ[it-2])×σk[:,ik,it-2]
                        -5.0/24.0*bk(k,Δ[it-1])×σk[:,ik,it-1]
                        +19.0/24.0*bk(k,Δ[it])×σk[:,ik,it]
                        +3.0/8.0*bk(k,Δ[it+1])×σk[:,ik,it+1])
                end
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
