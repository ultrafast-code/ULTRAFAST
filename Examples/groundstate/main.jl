#Main script to initialize:
#-number of workers for parallel computing
#-neural network and lattice
#-hyper-parameters of ground state and dynamic simulations

############################ Parallelizzation ##################################
using Distributed
#How many workers do you want? Set Nworkers to number of CPUs
#you want to add on top of the master CPU. Default Nworkers = 0
Nworkers = 3
addprocs(Nworkers)

################################################################################
###################### Neural network and lattice ##############################
#How large is your lattice? Set the lattice size Lx and Ly
Lx = 4
Ly = 4

#How many hidden units M=Alpha*Lx*Ly do you want? Set Alpha
Alpha = 2

#What is the exchange interaction along x and y? Set Jx and Jy (default = 1)
Jx = 1
Jy = 1

################################################################################
####################### Ground state optimization ##############################
#Do you want to run a ground state optimization? Set gs  (default = true)
gs = true

#Set hyperparameters of ground state optimization
#How many ground state iterations do you want? Set nIter
nIter = 300

#How many Monte Carlo samples do you want to use? Set nSampleGs
nSampleGs = 2000

#How small do you want the learning rate? Set stepSize
stepSize = 0.005

################################################################################
############################## Dynamics ########################################
#Do you want to run dynamics? Set dy (default dy = true)
dy = false

#What is the total integration time you want to simulate? Set totTime
totTime = 10.0

#How many Monte Carlo samples do you want to use? Set nSampleDyn
nSampleDyn = 1000

#How small do you want the step size of the Heun integrator? Set stepSizeHeun
stepSizeHeun = 0.001

#Do you want to invert the covariance matrix with MINRES or Pseudoinverse?
#Set minres_ = true to solve with MINRES, false to solve with Pseudoinverse (default = true)
minres_ = true

#What are the initial parameters for the dynamic W(t=0)? Import independent variational parameters
#Note that real and imaginary parts are conveniently stored in separate files
#First import package DelimitedFiles.jl to read and write files
using DelimitedFiles
#WRBMInit = readdlm("src/pre_optimized_parameters/W_RBM_4_2_real.jl") .+ im*readdlm("src/pre_optimized_parameters/W_RBM_4_2_imag.jl")

#What is the perturbation do you want to use along y? Define DeltaJ(t) (default DeltaJ(t)=0.1)
@everywhere function DeltaJ(t::Float64)
    return 0.1*exp(-(t - 0.987)^2 / (2*0.4^2))*sin(8*t)
end
################################# end ##########################################
################################################################################
include("../main/init.jl")
