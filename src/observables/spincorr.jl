#Calculate spin-spin correlation <S_1 \cdot S_2> for ground state wavefunction
function Spincorr_GS(WRBM::Array{Complex{Float64},1},Nsample::Int64,i_::Int64,j_::Int64)

    #Initialize network parameters
    RBM_Par = init_RBM_par(alpha,nspins)

    #Divide number of samples among CPU workers
    nsample = partition_sample(Nsample)

    #Store number of workers for later use
    NW = Int64[]; push!(NW, NumberWorkers())

    #Update neural network parameters at the current time t
    for j in procs()
        fetch(@spawnat j WToPar(RBM_Par,WRBM))
    end

    #Each worker evaluate spin-spin correlations and sum the outcome
    SP = @distributed (+) for y = 1:NW[1]
        sp_x(RBM_Par,y,nsample,i_,j_)
    end

    spincorr = SP./NW[1]

    return spincorr

end

#Calculate spin-spin correlation <S_1 \cdot S_2> for dynamic wavefunction
function Spincorr_DYN(WRBM_t::Array{Complex{Float64},2},Nsample::Int64,i_::Int64,j_::Int64)

    #Initialize network parameters
    RBM_Par = init_RBM_par(alpha,nspins)

    #Get time interval length from size of WRBM_t
    sizeWRBMt = length(WRBM_t[:,1])

    #Divide number of samples among CPU workers
    nsample = partition_sample(Nsample)

    #Store number of workers for later use
    NW = Int64[]; push!(NW, NumberWorkers())

    #Initialize array storing spin-spin correlations
    spincorr_t = zeros(sizeWRBMt)

    save_ind = 1
    for i=1:sizeWRBMt

        #Update neural network parameters at the current time t
        for j in procs()
            fetch(@spawnat j WToPar(RBM_Par,WRBM_t[i,:]))
        end

        #Each worker evaluate spin-spin correlations and sum the outcome
        SP = @distributed (+) for y = 1:NW[1]
            sp_x(RBM_Par,y,nsample,i_,j_)
        end

        spincorr_t[i] = SP./NW[1]

        #Every 10 steps correlations are saved for back-up
        save_ind += 1
        if save_ind == 10
            writedlm("Spincorr_$(nspins)_$(alpha).jl",spincorr_t)
            save_ind = 1
        end

    end

    return spincorr_t

end


#Evaluate <S_1 \cdot S_2> with Monte Carlo sampling
function sp_x(RBMPar::RBM_par,y::Int64,_nsample::Int64,i_::Int64,j_::Int64)

    #generate rabdom seed for the Monte Carlo chain depending on the worker y
    Random.seed!(100*y)

    #Set number of spin flips for the Monte Carlo chain. Default nflips=2
    nflips=SpinFlips()
    flips = Array{Int64}(undef,nflips)

    #Initialize a random state with zero magnetization if mag0=true
    state = GenRandomState(mag0)

    #initializing look-up tables in the wave-function
    Lt = InitLt(RBMPar,state)

    #thermalization
    @inbounds for n = 1:_nsample*0.1
      @inbounds for i = 1:nspins
                 StateUpdate(RBMPar,Lt,state,flips,nflips)
      end
    end

    spinc = Float64(0)
    #sequence of sweeps
    @inbounds for n = 1:_nsample
        @inbounds for i = 1:nspins
                    StateUpdate(RBMPar,Lt,state,flips,nflips)
                end
                spinc += real(sp_bond(RBMPar,Lt,state,[i_;j_]))
        end

        return spinc/_nsample

end
