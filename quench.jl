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
    σk = zeros(ComplexF64, 3, Nk, Nt+1);

    # Pauli matrices
    τ1 = [0 1; 1 0];
    τ2 = [0 -im; im 0];
    τ3 = [1 0; 0 -1];
    τμ = [τ1,τ2,τ3];
    
    # BdG Hamiltonian
    function hk(k_,Δ_)
        return -2*t_hop*cos(k_)*τ3-Δ_*(τ1-im*τ2)/2-conj(Δ_)*(τ1+im*τ2)/2;
    end

    # Effective magnetic field
    function bk(iμ_,k_,Δ_)
        if iμ_ == 1
            return -real(Δ_);
        elseif iμ_ == 2
            return -imag(Δ_);
        elseif iμ_ == 3
            return -2*t_hop*cos(k_);
        end
    end

    # Levi-Civita symbol
    function ϵ_ijk(iμ1_,iμ2_,iμ3_)
        if iμ1_ == 1 && iμ2_ == 2 && iμ3_ == 3
            return 1;
        elseif iμ1_ == 1 && iμ2_ == 3 && iμ3_ == 2
            return -1;
        elseif iμ1_ == 2 && iμ2_ == 1 && iμ3_ == 3
            return -1;
        elseif iμ1_ == 2 && iμ2_ == 3 && iμ3_ == 1
            return 1;
        elseif iμ1_ == 3 && iμ2_ == 1 && iμ3_ == 2
            return 1;
        elseif iμ1_ == 3 && iμ2_ == 2 && iμ3_ == 1
            return -1;
        else
            return 0;
        end
    end

    # equilibrium mean-field self-consistency
    diff = 1.0;
    iter = 0;
    while diff > tolerance
        Δ_temp = 0.0;
        for ik in 1:Nk
            k = ik*dk;
            Ek,Uk = eigen(hk(k,Δ[1]));
            f = diagm([
                1/(exp(Ek[1]/T)+1.0),
                1/(exp(Ek[2]/T)+1.0)]);
                Δ_temp += tr(f*Uk'*(τ1+im*τ2)/2*Uk);
        end
        Δ_temp = Δ_temp*Vi/Nk;

        diff = abs(Δ[1]-Δ_temp);
        Δ[1] = Δ_temp;
        iter += 1;
    end
    iter_list[1] = iter;
    Δ_re[1] = real(Δ[1]);
    Δ_im[1] = imag(Δ[1]);

    # initial pseudospin configuration
    E_tot[1] = 0.0;
    for ik in 1:Nk
        k = ik*dk;
        Ek,Uk = eigen(hk(k,Δ[1]));
        f = diagm([
            1/(exp(Ek[1]/T)+1.0),
            1/(exp(Ek[2]/T)+1.0)]);
        for iμ in 1:3
            σk[iμ,ik,1] = tr(f*Uk'*τμ[iμ]*Uk)/2;
        end
        for iμ in 1:3
            E_tot[1] += real(bk(iμ,k,Δ[1])*σk[iμ,ik,1]);
        end
    end
    E_tot[1] = 2*E_tot[1]/Nk+abs2(Δ[1])/Vi;

    # time evolution
    t_list[1] = 0.0;
    for it in 1:Nt
        t = it*dt;
        t_list[it+1] = t;
        Δ[it+1] = Δ[it];
        for ik in 1:Nk
            for iμ in 1:3
                σk[iμ,ik,it+1] = σk[iμ,ik,it];
            end
        end
    
        diff = 1.0;
        iter = 0;
        while diff > tolerance
            Δ_temp = 0.0;
            for ik in 1:Nk
                k = ik*dk;
                σk_temp = zeros(ComplexF64,3);
                for iμ1 in 1:3
                    σk_temp[iμ1] = σk[iμ1,ik,it];
                    for iμ3 in 1:3
                        for iμ2 in 1:3
                            if ϵ_ijk(iμ1,iμ2,iμ3) != 0
                                σk_temp[iμ1] += dt*ϵ_ijk(iμ1,iμ2,iμ3)*(
                                    bk(iμ2,k,Δ[it])*σk[iμ3,ik,it]
                                    +bk(iμ2,k,Δ[it+1])*σk[iμ3,ik,it+1]
                                );
                            end
                        end
                    end
                end
                Δ_temp += σk_temp[1]+im*σk_temp[2];
                for iμ in 1:3
                    σk[iμ,ik,it+1] = σk_temp[iμ];
                end
            end
            Δ_temp = Δ_temp*Vf/Nk;
            diff = abs(Δ[it+1]-Δ_temp);
            Δ[it+1] = Δ_temp;
            iter += 1;
        end
        iter_list[it+1] = iter;
        Δ_re[it+1] = real(Δ[it+1]);
        Δ_im[it+1] = imag(Δ[it+1]);
        E_tot[it+1] = 0.0;
        for ik in 1:Nk
            k = ik*dk;
            for iμ in 1:3
                E_tot[it+1] += real(bk(iμ,k,Δ[it+1])*σk[iμ,ik,it+1]);
            end
        end
        E_tot[it+1] = 2*E_tot[it+1]/Nk+abs2(Δ[it+1])/Vf;
    end
end

time_dep_BCS()

fig_Delta_re = plot(t_list,[Δ_re,Δ_im],xlim=(0,50),ylim=(0,4),xlabel="t",ylabel="Δ",label=["Re Δ" "Im Δ"])

fig_E_tot = plot(t_list,E_tot,xlim=(0,50),ylim=(-1.5,-1.3),xlabel="t",ylabel="E",label="E_tot")

savefig(fig_Delta_re,"Delta_re.png")

savefig(fig_E_tot,"E_tot.png")
