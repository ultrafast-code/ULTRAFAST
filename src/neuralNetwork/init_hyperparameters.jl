#initialize hyper-parameters of ground state and dynamic optimizations

#define structure type containing hyper-parameters of ground state optimization:
#number of Monte Carlo samples 'samples'
#number of iteration steps 'n_iter'
#size of the learning rate 'step_size'
struct GS_hp
    samples::Int64
    n_iter::Int64
    step_size::Float64
end

#define structure type containing hyper-parameters of dynamic optimization:
#number of Monte Carlo samples 'samples'
#total integration time 'tot_time'
#step size of the time-integration scheme 'step_size'
struct DYN_hp
    samples::Int64
    tot_time::Float64
    step_size::Float64
    minres_::Bool
end
