using LinearAlgebra, Plots

Nt = 1000;
t_list = zeros(Float64, Nt+1);
Δ_re = zeros(Float64, Nt+1);
Δ_im = zeros(Float64, Nt+1);
E_tot = zeros(Float64, Nt+1);
iter_list = zeros(Int16, Nt+1);

function time_dep_BCS()
    # parameters
    Vi = 2.0;
    Vf = 6.0;
    t_hop = 1.0;
    T = 0.05;
    Nk = 200;
    dk = 2.0*π/Nk;
    tmax = 50.0;
    dt = tmax/Nt;
    tolerance = 0.000001;

    # gap function
    Δ = zeros(ComplexF64, Nt+1);
    Δ[1] = 1.0;

    # pseudospin
    σk = zeros(Float64, 3, Nk, Nt+1);

    # Pauli matrices
    τ1 = [0 1; 1 0];
    τ2 = [0 -im; im 0];
    τ3 = [1 0; 0 -1];
    τ = [τ1,τ2,τ3];
    
    # BdG Hamiltonian
    function hk(k::Float64, Δ::ComplexF64)
        return -2*t_hop*cos(k)*τ3-Δ*(τ1-im*τ2)/2-conj(Δ)*(τ1+im*τ2)/2;
    end

    # Effective magnetic field
    function bk(k::Float64, Δ::ComplexF64)
        return [-real(Δ),-imag(Δ),-2*t_hop*cos(k)];
    end

    # equilibrium mean-field self-consistency
    diff = 1.0;
    iter = 0;
    while diff > tolerance
        Δ_temp = 0.0;
        for ik in 1:Nk
            k = ik*dk;
            Ek,Uk = eigen(hk(k,Δ[1]));
            f = diagm([1/(exp(Ek[1]/T)+1.0), 1/(exp(Ek[2]/T)+1.0)]);
            Δ_temp += tr(f*Uk'*(τ1+im*τ2)/2*Uk);
        end
        Δ_temp = Δ_temp*Vi/Nk;
        diff = abs(Δ[1]-Δ_temp);
        Δ[1] = Δ_temp;
        iter += 1;
    end
    t_list[1] = 0.0;
    Δ_re[1] = real(Δ[1]);
    Δ_im[1] = imag(Δ[1]);
    iter_list[1] = iter;

    # initial pseudospin configuration
    E_tot[1] = 0.0;
    for ik in 1:Nk
        k = ik*dk;
        Ek,Uk = eigen(hk(k,Δ[1]));
        f = diagm([1/(exp(Ek[1]/T)+1.0), 1/(exp(Ek[2]/T)+1.0)]);
        for α in 1:3
            σk[α,ik,1] = real(tr(f*Uk'*τ[α]*Uk)/2);
        end
        E_tot[1] += real(bk(k,Δ[1])⋅σk[:,ik,1]);
    end
    E_tot[1] = 2*E_tot[1]/Nk + abs2(Δ[1])/Vi;

    # time evolution
    for it in 1:Nt
        t = it*dt;
        Δ[it+1] = Δ[it];
        σk[:,:,it+1] = σk[:,:,it];
        diff = 1.0;
        iter = 0;
        while diff > tolerance
            Δ_temp = 0.0;
            for ik in 1:Nk
                k = ik*dk;
                σk[:,ik,it+1] = σk[:,ik,it] + dt*(bk(k,Δ[it])×σk[:,ik,it]+bk(k,Δ[it+1])×σk[:,ik,it+1]);
                Δ_temp += σk[1,ik,it+1] + im*σk[2,ik,it+1];
            end
            Δ_temp = Δ_temp*Vf/Nk;
            diff = abs(Δ[it+1]-Δ_temp);
            Δ[it+1] = Δ_temp;
            iter += 1;
        end
        t_list[it+1] = t;
        Δ_re[it+1] = real(Δ[it+1]);
        Δ_im[it+1] = imag(Δ[it+1]);
        iter_list[it+1] = iter;
        E_tot[it+1] = 0.0;
        for ik in 1:Nk
            k = ik*dk;
            E_tot[it+1] += real(bk(k,Δ[it+1])⋅σk[:,ik,it+1]);
        end
        E_tot[it+1] = 2*E_tot[it+1]/Nk + abs2(Δ[it+1])/Vf;
    end
end

time_dep_BCS()

fig_Delta = plot(t_list,[Δ_re,Δ_im],xlim=(0,50),ylim=(0,4),xlabel="t",ylabel="Δ",label=["Re Δ" "Im Δ"])
fig_E_tot = plot(t_list,E_tot,xlim=(0,50),ylim=(-1.5,-1.3),xlabel="t",ylabel="E",label="E_tot")

savefig(fig_Delta,"Delta.png")
savefig(fig_E_tot,"E_tot.png")