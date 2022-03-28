using IterativeSolvers
dir = Base.source_dir()
@everywhere include(joinpath($dir,"model.jl"))
################################################################################
@everywhere include(joinpath($dir,"initialization.jl"))

########## GROUND STATE OPTIMIZATION ##############################
nsweeps = 1000                   #number of monte carlo samples
n_iter = 200                     #number of iterations
gamma = 0.005                    #step size during optimization

gs_optimization(n_iter,nsweeps,gamma)                 #run a ground state optimization
writedlm("W_rbm_$(nspins)_$(nhv[1]).jl",W_RBM)        #save W_RBM in "W_rbm_nspins_alpha.jl"

nsweeps_gs = 10000                           #number of samples for g.s. observables evaluation
m = 1; n = 2;                                #select indices of spin-spin correlation function <S_m S_n>
GS_obs(W_RBM,nsweeps_gs,m,n)                 #calculate energy per spin and <S m S n > in the g.s.
########### UNITARY DYNAMICS ###################################################
nsweeps_d = 2000                             #number of Monte Carlo sweeps in the time-evolution
length_int = 1.                              #time-length of the time integrationo
step_size = 0.0025                           #time-step of the time integration

Init = W_RBM                                 #set the optimized parameters as initial wavefunction
run_dynamics(heun,Init)                      #run the time-evolution with initial parameters Init
############ OBSERVABLE EVALUATION #############################################
nsweeps_obs = 10000                          #number of samples for the evaluation of <S i S j>
i = 1; j = 2;                                #select indices of spin-spin correlation function <S_i S_j>

Spincorr_d = spincorr_d(W_RBM_t,nsweeps_obs,step_size,i,j);;      #evaluate the time-evolution of  <S_i S_j>
