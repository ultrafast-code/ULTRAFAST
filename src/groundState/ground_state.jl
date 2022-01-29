#Define main function for the ground state optimization
#and import the necessary routines

#Main function for ground state optimization
function GS_optimization(GsHp::GS_hp, dim::Int64)

    println("# Starting ground state optimization for the Heisenberg Hamiltonian, with $nspins = $L_x x $L_y spins and α=$alpha")
    println("# Number of workers = $(nprocs())")
    println("# The following hyper-parameters are used:")
    println("# Number of sweeps = $(GsHp.samples),  iteration step = $(GsHp.n_iter),  learning rate = $(GsHp.step_size)")


    #If the number of workers exceeds the number of variational parameters,
    #only a number of workers equals to the number of variational parameters will be used.
    #The exceeding workers will be removed.
    if nprocs() > dim
        println("WARNING: number of workers exceeds number of variational parameters. Only $dim workers will be used.")
        diff = nprocs() - dim
        rmprocs(nprocs()-diff+1:nprocs())
    end

    #Store number of workers for later use
    NW = Int64[]; push!(NW, NumberWorkers())

    #Define matrix storing energy and wavefunction gradients for each sampled state
    #@everywhere ensure a copy of ENWfDer is created for each worker
    @everywhere ENWfDer = Array{Complex{Float64}}(undef,@getfrom 1 Nsample,dim + 1)

    #Initialize network parameters to random numbers
    RBM_Par = init_RBM_par(alpha,nspins)
    #Sort network parameters to W_RBM independent parameters and vice versa
    W_RBM = Array{Complex{Float64}}(undef,dim)
    ParToW(RBM_Par,W_RBM); WToPar(RBM_Par,W_RBM)

    #Set value of J_x, J_y
    @everywhere J = [J_x; J_y]
################################################################################
#############    Start main iteration  #########################################
    for NIter = 1:GsHp.n_iter

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

         #Regularization of the S-matrix with a small shift of the diagonal element of S
         SMatrix = SMatrix .+ Matrix{Float64}(I,dim,dim).*max(100*(0.9)^NIter,0.0001)

         #Inversion of the S-matrix with normal inversion routine (not pseudo-inverse)
         InvS = inv(SMatrix)

         #Store the total energy for each state
         energy = EnergyStates[:]
         #Calculate the mean value of the energy and energy variance and store them
         mean_E = mean(energy); var_E = var(energy)
         push!(Energy_,real(mean_E)); push!(Variance_,norm(var_E))

         ##calculate the derivative of the energy respect to the variational parameters and assign them to the array dF
         EGrad = @distributed (+) for i=1:NW[1]
            Egrad(dim)./NW[1]
         end

         #Update variational parameters according to SR optimization rule
         for i=1:dim
             W_RBM[i] = W_RBM[i] - GsHp.step_size*(transpose(InvS[i,:])*EGrad)
         end

         #Update neural network parameters after optimization
         for i in procs()
             fetch(@spawnat i WToPar(RBM_Par,W_RBM))
         end

         #print value of energy every 10 steps
         if NIter % 10 == false
             println("Iteration step #$NIter")
             println("Energy: $(real(mean_E)) +- $(norm(var_E))")
             println(" ")
         end

     end

###########################   End main iteration    ############################
################################################################################
    println("# ground state optimization completed")
    println(" ")

    #Store the final independent variational parameters (real and imaginary parts)
    writedlm("W_RBM_$(nspins)_$(alpha)_real.jl", real(W_RBM)); writedlm("W_RBM_$(nspins)_$(alpha)_imag.jl", imag(W_RBM))

end
