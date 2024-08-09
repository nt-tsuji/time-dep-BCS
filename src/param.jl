# parameters
struct Param
    Nt::Int64           # number of time steps
    tmax::Float64       # maximum time
    Vi::Float64         # initial interaction strength
    Vf::Float64         # final interaction strength
    t_hop::Float64      # hopping amplitude
    T::Float64          # initial temperature
    Nk::Int64           # number of k points
    tolerance::Float64  # tolerance of convergence
    dk::Float64         # dk = 2π/Nk
    dt::Float64         # dt = tmax/Nt
    Param(Nt, tmax, Vi, Vf, t_hop, T, Nk, tolerance) = new(Nt, tmax, Vi, Vf, t_hop, T, Nk, tolerance, 2.0*π/Nk, tmax/Nt)
end