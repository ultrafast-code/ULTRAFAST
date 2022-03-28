const pert = Array{Bool}(1)
const pert[1] = false
typer = 0*im
const coupling= Array{Float64}(1)
const coupling[1]=0.1

if isa(α,Int64)==false || α <= 0
    error("α must be a positive integer")
end

if mag0 == true && n_flips == 1
    error("Zero magnetization sector not compatible with n_flips=1")
end

const nhv = Array{Int64}(2)
const nhv[:] = [α,n_spins]
const heis_2D = false
const jz = 1.
Ham = "HH_2D"

if round(sqrt(nhv[2]))- (sqrt(nhv[2])) != 0. || rem(nhv[2],2) != 0
    error("The number of spins must be N=LxL with L even.")
end
include("init_functions.jl")
include("2D_heisenberg.jl")



 if Sym == "No symmetry"
  include("init_nosymmetry.jl")
elseif Sym == "Translation symmetry" && pbc == true
  include("init_translinv.jl")
 else error("select symmetry compatible with the chosen Hamiltonian.")
end
InitVarPar(nspins,1049)
include("observables.jl")
