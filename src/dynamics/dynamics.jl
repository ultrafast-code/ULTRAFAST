
#Start real-time evolution from initial variational parameters WRBM_Init and with "method" integration scheme
function runDynamics(method,DYNHp::DYN_hp,WRBM_Init::Array{Complex{Float64},1})

    println("# Starting dynamics for $nspins = $L_x x $L_y spins and α=$alpha")
    println("# Time evolution = [0.,$(DYNHp.tot_time)], time-step = $(DYNHp.step_size). ")
    println("# Number of sweeps = $(DYNHp.samples)")
    println("# Number of workers = $(nprocs())")

    #Define array storing independent variational parameters at each time step
    #starting from the initial parameters WRBM_Init
    WRBM_t = Array{Array{Complex{Float64}}}(undef,0)
    push!(WRBM_t, WRBM_Init)

    #Iniziate the unitary time evolution
    @time method(RHS,WRBM_t,DYNHp)

end

#Heun integration scheme
function heun(GenRHS,WRBMt::Array{Array{Complex{Float64}}},DYNHp::DYN_hp)

    #initialized the time-interval of the simulation
    time = collect(0.:DYNHp.step_size:DYNHp.tot_time)

    #Assign independent variational parameters at t=0 to new variable W_
    W_ = WRBMt[1][:]

    #Initialize network parameters
    RBM_Par = init_RBM_par(alpha,nspins)
    #Sort network parameters to W_RBM independent parameters
    WToPar(RBM_Par,W_)

    #If the number of workers exceeds the number of variational parameters,
    #only a number of workers equals to the number of variational parameters will be used.
    #The exceeding workers will be removed.
    if nprocs() > dim
        println("WARNING: number of workers exceeds number of viariational parameters. Only $dim workers will be used.")
        diff = nprocs()-dim
        rmprocs(nprocs()-diff+1:nprocs())
    end

    #Define matrix storing energy and wavefunction gradients for each sampled state
    #@everywhere ensure a copy of ENWfDer is created for each worker
    @everywhere ENWfDer = Array{Complex{Float64}}(undef,@getfrom 1 NsampleDYN,dim + 1)

    save_ind = 1
    for t in time
            RHS_t = GenRHS(WRBMt[end],t,RBM_Par,true,DYN_HP.minres_)
            W_ .+=  DYNHp.step_size*(RHS_t .+ GenRHS(WRBMt[end] .+ DYNHp.step_size*RHS_t,t+DYNHp.step_size,RBM_Par,false,DYN_HP.minres_))/2
            push!(WRBMt, W_[:] ); save_ind += 1

            #Backup save of independent variational parameters every ten time steps
            if save_ind == 10
                writedlm("W_RBM_t_$(nspins)_$(alpha)_real.jl",real(WRBMt));writedlm("W_RBM_t_$(nspins)_$(alpha)_imag.jl",imag(WRBMt))
                save_ind = 1
            end

    end

    #Save indipendent variational parameters at each time step t
    writedlm("W_RBM_t_$(nspins)_$(alpha)_real.jl",real(WRBMt));writedlm("W_RBM_t_$(nspins)_$(alpha)_imag.jl",imag(WRBMt))

end

#Generate the right hand side of the TDVP equations  d/dt W(t) = -i*S^-1 F(W(t))
function RHS(W_t::Array{Complex{Float64},1},t::Float64,RBM_Par::RBM_par,time_flag::Bool,minres_::Bool)

    #Store number of workers for later use
    NW = Int64[]; push!(NW, NumberWorkers())

    #Define tt to pass current time to each worker
    @everywhere tt = Array{Float64}(undef,1)
    for i in procs()
    fetch(@spawnat i tt[1] = t )
    end

    #perturbation of J_ex along y
    @everywhere J = [J_x + DeltaJ(tt[1]); J_y]

    #Update neural network parameters at the current time t
    for i in procs()
        fetch(@spawnat i WToPar(RBM_Par,W_t))
    end

    #Sampling of states
    #At the end of the loop, the energy sampled are returned in EnergyStates
    #The derivatives of the wavefunction are stored in ENWFDer
     EnergyStates = @distributed (vcat) for i=1:NW[1]
          sampling(RBM_Par,J,i)
     end

     #Evaluation entries of S-matrix stored in the vector SVector
     SVector = @distributed (+) for i=1:NW[1]
        S_gen(dim)./NW[1]
     end

     #Sort entries SVector into the S matrix
     SMatrix = Gen_S(SVector)

     #Store the total energy for each state
     energy = EnergyStates[:]
     #Calculate the mean value of the energy and energy variance
     mean_E = mean(energy); var_E = var(energy)

     #Calculate the derivative of the energy respect to the variational parameters and assign them to the array dF
     EGrad = @distributed (+) for i=1:NW[1]
        Egrad(dim)./NW[1]
     end

     #Print the energy only in the first step of the Heun scheme
      if time_flag == true
      println("time $t- total energy: $(real(mean_E)) +- $(real(var_E)/length(energy))")
      end

      #Return the rhs of the TDVP equation
      if minres_ == true
          return rhsMinres(SMatrix,EGrad)
      else
          return rhsPinv(SMatrix,EGrad)
      end

end

#Generate rhs of TDVP equation by using the iterative solver MINRES
function rhsMinres(S_::Matrix{Complex{Float64}}, eGrad::Array{Complex{Float64},1})

    return minres(S_,-im*eGrad)

end

#Generate rhs of TDVP equation by using singular value decomposition (pseudo inverse)
function rhsPinv(S_::Matrix{Complex{Float64}}, eGrad::Array{Complex{Float64},1})

    #Regularization of the S-matrix with a small shift of the diagonal element of S
    S_ .= S_ .+ Matrix{Float64}(I,dim,dim).*0.0001

    #Calculate the pseudo inverse of S_
    invS = pinv(S_)

    return -im*invS*eGrad

end


#Tamed Euler integration scheme
function tamedEuler(GenRHS,WRBMt::Array{Array{Complex{Float64}}},DYNHp::DYN_hp)

    #initialized the time-interval of the simulation
    time = collect(0.:DYNHp.step_size:DYNHp.tot_time)

    #Assign independent variational parameters at t=0 to new variable W_
    W_ = WRBMt[1][:]

    #Initialize network parameters
    RBM_Par = init_RBM_par(alpha,nspins)
    #Sort network parameters to W_RBM independent parameters
    WToPar(RBM_Par,W_)

    #If the number of workers exceeds the number of variational parameters,
    #only a number of workers equals to the number of variational parameters will be used.
    #The exceeding workers will be removed.
    if nprocs() > dim
        println("WARNING: number of workers exceeds number of viariational parameters. Only $dim workers will be used.")
        diff = nprocs()-dim
        rmprocs(nprocs()-diff+1:nprocs())
    end

    #Define matrix storing energy and wavefunction gradients for each sampled state
    #@everywhere ensure a copy of ENWfDer is created for each worker
    @everywhere ENWfDer = Array{Complex{Float64}}(undef,@getfrom 1 NsampleDYN,dim + 1)

    save_ind = 1
    for t in time
            RHS_t = GenRHS(WRBMt[end],t,RBM_Par,true,DYNHp.minres_)
            W_ .+= (RHS_t * DYNHp.step_size)./(1 .+ DYNHp.step_size*norm(RHS_t))
            push!(WRBMt, W_[:] ); save_ind += 1

            #Backup save of independent variational parameters every ten time steps
            if save_ind == 10
                writedlm("W_RBM_t_$(nspins)_$(alpha)_real.jl",real(WRBMt));writedlm("W_RBM_t_$(nspins)_$(alpha)_imag.jl",imag(WRBMt))
                save_ind = 1
            end

    end

    #Save independent variational parameters at each time step t
    writedlm("W_RBM_t_$(nspins)_$(alpha)_real.jl",real(WRBMt));writedlm("W_RBM_t_$(nspins)_$(alpha)_imag.jl",imag(WRBMt))

end

#Tamed Heun integration scheme
function tamedHeun(GenRHS,WRBMt::Array{Array{Complex{Float64}}},DYNHp::DYN_hp)

    #initialized the time-interval of the simulation
    time = collect(0.:DYNHp.step_size:DYNHp.tot_time)

    #Assign independent variational parameters at t=0 to new variable W_
    W_ = WRBMt[1][:]

    #Initialize network parameters
    RBM_Par = init_RBM_par(alpha,nspins)
    #Sort network parameters to W_RBM independent parameters
    WToPar(RBM_Par,W_)

    #If the number of workers exceeds the number of variational parameters,
    #only a number of workers equals to the number of variational parameters will be used.
    #The exceeding workers will be removed.
    if nprocs() > dim
        println("WARNING: number of workers exceeds number of viariational parameters. Only $dim workers will be used.")
        diff = nprocs()-dim
        rmprocs(nprocs()-diff+1:nprocs())
    end

    #Define matrix storing energy and wavefunction gradients for each sampled state
    #@everywhere ensure a copy of ENWfDer is created for each worker
    @everywhere ENWfDer = Array{Complex{Float64}}(undef,@getfrom 1 NsampleDYN,dim + 1)

    save_ind = 1
    for t in time
            RHS_t1 = GenRHS(WRBMt[end],t,RBM_Par,true,DYNHp.minres_)
            RHS_t2 = GenRHS(WRBMt[end] .+ DYNHp.step_size*RHS_t1,t+DYNHp.step_size,RBM_Par,false,DYN_HP.minres_)
            W_ .+=  0.5*DYNHp.step_size*(RHS_t1/(1 + DYNHp.step_size*norm(RHS_t1) ) .+  RHS_t2/(1 + DYNHp.step_size*norm(RHS_t2)))
            push!(WRBMt, W_[:] ); save_ind += 1

            #Backup save of independent variational parameters every ten time steps
            if save_ind == 10
                writedlm("W_RBM_t_$(nspins)_$(alpha)_real.jl",real(WRBMt));writedlm("W_RBM_t_$(nspins)_$(alpha)_imag.jl",imag(WRBMt))
                save_ind = 1
            end

    end

    #Save independent variational parameters at each time step t
    writedlm("W_RBM_t_$(nspins)_$(alpha)_real.jl",real(WRBMt));writedlm("W_RBM_t_$(nspins)_$(alpha)_imag.jl",imag(WRBMt))

end
