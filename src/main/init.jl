#Initiailize ULTRAFAST according to the parameters set in main.jl
#Note that the use of const global variable is more efficient than global variables

@everywhere using ParallelDataTransfer; @everywhere using DelimitedFiles; @everywhere using IterativeSolvers
@everywhere using Statistics; @everywhere using LinearAlgebra
#Define as const global variables the parameters set in main.jl
@everywhere const alpha = Int64(@getfrom 1 Alpha)
@everywhere const L_x = Int64(@getfrom 1 Lx)
@everywhere const L_y = Int64(@getfrom 1 Ly)
@everywhere const nspins = L_x*L_y            #Number of spins in the lattice

@everywhere const J_x = Float64(@getfrom 1 Jx)
@everywhere const J_y = Float64(@getfrom 1 Jy)
@everywhere const dim = alpha*nspins + alpha  #Number of independent variational parameters

#Set magnetization sector for Monte Carlo sampling. (default mag0 = true, only sample zero magnetization states)
@everywhere const mag0 = true

@everywhere include("../neuralNetwork/init_hyperparameters.jl")
@everywhere include("../neuralNetwork/nn.jl")
@everywhere include("../neuralNetwork/WtoPar_ti.jl")
@everywhere include("../observables/S_matrix.jl")
@everywhere include("../observables/GradE.jl")
@everywhere include("../parallelization/parallelization.jl")
@everywhere include("../observables/energy.jl")
@everywhere include("../observables/wv_derivative.jl")
@everywhere include("../sampling/sampling.jl")
@everywhere include("../lattice/Heisenberg_lattice.jl")
@everywhere include("../groundState/ground_state.jl")
@everywhere include("../dynamics/dynamics.jl")
@everywhere include("../observables/spincorr.jl")
@everywhere cd("..\\ULTRAFAST\\src\\output")
#Define struct storing hyper-parameters set in main.jl
GS_HP = GS_hp(nSampleGs,nIter,stepSize)
DYN_HP = DYN_hp(nSampleDyn,totTime,stepSizeHeun,minres_)

#Store nearest neighbouring bonds
@everywhere const bonds_x, bonds_y = Lattice.InitLattice(L_x,L_y,nspins)

#Store filter matrix for efficient evaluation of wavefunction derivatives
@everywhere const Filter_ind = WvDerivativeTi.transl_inv_tool(nspins,alpha,L_x,L_y)

#Define number of samples handled by each worker
@everywhere const Nsample = partition_sample(@getfrom 1 GS_HP.samples)
@everywhere const NsampleDYN = partition_sample(@getfrom 1 DYN_HP.samples)

############################### Start simulations ##############################

#Start ground state optimization
if gs == true
    const Energy_ = Float64[]; const Variance_ = Float64[]
    @time GS_optimization(GS_HP,dim)
end

#Start dynamics
if dy == true
    runDynamics(tamedHeun,DYN_HP,WRBMInit[:])
end
