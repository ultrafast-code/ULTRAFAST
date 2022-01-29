#Define functions for the Monte Carlo sampling:: Symmetry independent

#Sample nsweeps states and calculate energy and wavefunction derivatives for all the states
function sampling(RBMPar::RBM_par,_J::Array{Float64,1},y::Int64;thermfactor=0.1,sweepfactor=1)
       deriv_ = ENWfDer
       #Define number of samples from length of deriv_
       nsample = length(deriv_[:,1])

       #generate rabdom seed for the Monte Carlo chain depending on the worker y
       Random.seed!(100*y)

       #Set number of spin flips for the Monte Carlo chain. Default nflips=2
       nflips=SpinFlips()
       flips = Array{Int64}(undef,nflips)

       #Define quantity storing energy
       energy = Array{Complex{Float64}}(undef,0)

       #Initialize a random state with zero magnetization if mag0=true
       state = GenRandomState(mag0)

       #initializing look-up tables in the wave-function
       Lt = InitLt(RBMPar,state)

       #thermalization
       @inbounds for n=1:(nsample)*thermfactor
         @inbounds for i=1:(nspins*sweepfactor)
                    StateUpdate(RBMPar,Lt,state,flips,nflips)
         end
       end

       #sequence of sweeps to calculate energy and gradient of wavefunction
       @inbounds for n=1:nsample
           @inbounds    for i=1:(nspins*sweepfactor)
               StateUpdate(RBMPar,Lt,state,flips,nflips)
           end;

           Energy(RBMPar,Lt,state,energy,_J)
           WvDerivativeTi.Var_derivatives(state,deriv_,Lt,Filter_ind,nspins,alpha,n)

       end;

       @views deriv_[:,1] = energy

       return energy

end

#Initializes a random spin state with nspins size
#if mag0=true, the initial state is prepared with zero total magnetization
function GenRandomState(mag0::Bool)

      state_ = Int64[]

      for i=1:nspins
        push!(state_, 2*rand(0:1)-1)
      end;

      if mag0
        magt=1;
        if nspins%2 == 1
          print("# Error : Cannot initializate a random state with zero magnetization for odd number of spins")
          error()
        end;
        while magt!=0
          magt = 0;
          for i=1:nspins
            magt += state_[i]
          end;
          if magt>0
            rs = rand(1:nspins)
            while state_[rs]<0
              rs = rand(1:nspins)
            end;
            state_[rs] = -1
            magt -= 1

          elseif magt<0
             rs = rand(1:nspins)
            while state_[rs]>0
              rs = rand(1:nspins)
            end;
            state_[rs] = 1
            magt += 1

          end
        end
      end
      return state_
end


#Select the number of spin flips
function SpinFlips()
    return 2
end


#initialization of the look-up tables
function InitLt(RBMPar::RBM_par,state::Array{Int64,1})

  Lt_ = Array{Complex{Float64}}(undef,alpha*nspins)

  for h=1:alpha*nspins
    Lt_[h] = RBMPar.b[h]
    for v=1:nspins
      Lt_[h] += state[v]*RBMPar.W[v,h]
    end
  end

  return Lt_

end

#Generate a spin-flip in the given state
#the vector "flips" contains the indices of sites to be flipped
function GenSpin(state::Array{Int64,1},flips::Array{Int64,1},nflips::Int64,mag0)

      resize!(flips,nflips)

      flips[1]= rand(1:nspins)
      if nflips===2
        flips[2]=rand(1:nspins)
        if !mag0
          return flips[2]!== flips[1]
        else
          return state[flips[2]]!== state[flips[1]]
        end
      end

      return true

end


#Generate a new state from the given stat through a flip of the spins
#according to the acceptance probability and then update look-up table
#the vector "flips" contains the indices of sites to be flipped
function StateUpdate(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1},flips::Array{Int64,1},nflips::Int64)

    #Picking "nflips" random spins to be flipped
    if GenSpin(state,flips,nflips,mag0)

      #Computing acceptance probability
      acceptance = abs2(ExpLog(RBMPar,Lt,state,flips))

      #Metropolis-Hastings test
      if acceptance>rand()

        #Updating look-up tables in the wave-function
        UpdateLt(RBMPar,Lt,state,flips)

        #Moving to the new configuration
        for i=1:length(flips)
          state[flips[i]] *= -1
        end

      end
    end

end

#updates the look-up tables after spin flips
#the vector "flips" contains the indices of sites to be flipped
function UpdateLt(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1},flips::Array{Int64,1})

  if length(flips)===0
    return
  end

  for h=1:alpha*nspins
    for i=1:length(flips)
      Lt[h] -= 2.0*state[flips[i]]*RBMPar.W[flips[i],h]
    end
  end

end


#computes the logarithm of Psi(state')/Psi(state) where
#state' is a state with a certain number of flipped spins
# the vector "flips" contains the sites to be flipped
function Log_wf(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1},flips::Array{Int64,1})

      if length(flips) === 0
          return Complex(0)
      end;

      logwf = Complex(0)

      #Change due to the visible bias
        for i=1:length(flips)
          logwf -= RBMPar.a[flips[i]]*2*state[flips[i]]
        end;

        #Change due to the weights
@fastmath for  h=1:alpha*nspins
           thetah = Lt[h]
           thetahp = thetah

           for i=1:length(flips)
               thetahp -= 2.0*state[flips[i]]*RBMPar.W[flips[i],h]
           end

           logwf += log(cosh(thetahp))-log(cosh(thetah))
         end

         return logwf
end

#Return the exp of Log_wf()
function ExpLog(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1},flips::Array{Int64,1})
    return exp(Log_wf(RBMPar,Lt,state,flips))
end
