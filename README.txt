###################### Copyright information ##################################
ULTRAFAST is provided by G. Fabiani and J.H. Mentink as Supplementary Material
to the published paper:

"Investigating ultrafast quantum magnetism with machine learning"

ULTRAFAST is an open-source code for the simulation of ground state and dynamic
properties of quantum spin systems. Permission for using and coping part of or
all this code is granted, provided that a full citation of the paper published
with ULTRAFAST is displayed.

ULTRAFAST is written for research purposes only and distributed to reproduce the
results obtained in the paper. The use of this code is at the user's own risk.

######################  Running ULTRAFAST   ##################################
To run ULTRAFAST, first install Julia, version 0.6.1. Then install the code by 
downloading the tar file given in the Supplementary Material and untarring it in 
a suitable working directory. Julia can be run either from an interactive session
Read-Eval-Print Loop ("REPL")  or from the command line. To execute ULTRAFAST
in the REPL, double-click the Julia executable and type

julia> include("run.jl")

From the command line, open a terminal and type

$ path/to/julia   path/to/run.jl

The code features parallel computation. To run ULTRAFAST over N processes on the REPL, type

julia> addprocs(N)
julia> include("run.jl")

while on the command line, simply type

$ path/to/julia  -pN  path/to/run.jl

To start a simulation, the neural network and the physical problem to solve need
to be initialized. This can be done in the file "model.jl".
In the script "run.jl" you can choose to run both the ground state and the
dynamic optimization:
(1) Ground state optimization starts by calling the function
    gs_optimization(), which requires the number of Monte Carlo samples, the number
    of iterations and the learning rate as input. At the end of the optimization the
    optimal parameters W are stored in the file "W_rbm_nspins_alpha.jl".

(2) Dynamics is run by calling run_dynamics(). Analogous to the ground state
    optimization, this requires the number of Monte Carlo samples, the total
    evolution time and the time-step of the numerical time-integration as input.
    The variable "Init" is an array with the parameters of the initial state
    wavefunction. It can be initialized either by using the optimal parameters found
    in the ground state optimization (default option) or by choosing one of the
    pre-optimized wavefunction given in the Supplementary Information, which can be used
    only for translation invariant wave functions. For the latter define: 
    Init = readdlm("W_RBM_nspins_alpha_ti.jl"), where nspins (number of spins) 
    and alpha must be chosen according to "model.jl".

The functions GS_obs() and and spincorr_d() allow to measure < S_i S_j > for any
given i and j respectively in the ground state or along the time-evolution.
The function GS_obs() also provides the ground state energy per spin.

Further information about using the code can be found in the paper.