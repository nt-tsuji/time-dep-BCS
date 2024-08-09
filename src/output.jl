# output data
mutable struct Output
    t_list::Array{Float64,1}   # list of time steps
    Δ::Array{ComplexF64,1}     # gap function
    E_tot::Array{Float64,1}    # total energy
    iter_list::Array{Int16,1}  # list of numbers of iterations
    Output(Nt::Int64) = new(zeros(Float64,Nt+1), zeros(ComplexF64,Nt+1), zeros(Float64,Nt+1), zeros(Int64,Nt+1))
end