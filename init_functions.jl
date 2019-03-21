#This script contains immutable functions
#Define variational parameters
  const typer_bool = true
  const a = Array{Complex{Float64}}(40)   #visible bias
  const b = Array{Complex{Float64}}(40)   #hidden bias
  const W = Array{Array{Complex{Float64}}}(40) #weigths
  const Lt = Array{Complex{Float64}}(40)  #Look-up table
  W_RBM = Array{Complex{Float64}}(0)      #Array containing the actual variational parameters to be optimized

  const M = nhv[1]*nhv[2]                   #number of weigths
  const nspins = nhv[2]                     #number of spins
  const L=Int(sqrt(nhv[2]))                 #lattice size
  const l=Int(L/2)                          #half of the lattice size
  const flips = Array{Int64}(40);           #Array storing spin flips
  const W_RBM_t = Array{Array{Complex{Float64}}}(0)       #store the variational parameters at each time-step
  state = Array{Int64}(40);                 #Array storing the configuration state

#Main function for the ground state optimization
function gs_optimization(n_iterations::Int64,nsweeps::Int64,gamma::Float64,thermfactor=0.1,sweepfactor=1)
      println("# Starting optimization for $Ham, nspins=$nspins and α=$α. $Sym")
      println("# Number of workers = $(nprocs())")
      println("# Number of sweeps = $nsweeps")

      nflips=SpinFlips();

    #checking input consistency
    if nflips>2 || nflips<1
      print("# Error : The number of spin flips should be equal to 1 or 2. \n")
      error();
    end;
    if thermfactor>1 || thermfactor<0
      print("# Error : The thermalization factor should be a real number between 0 and 1 \n");
      error();
    end;
    if nsweeps<1
      print("# Error : Please enter a number of sweeps sufficiently large (>50)\n");
      error();
    end;

    #Initialize a random state with zero magnetization if mag0=true
    GenRandomState(mag0);

    resize!(flips,nflips);

    #initializing look-up tables in the wave-function
    InitLt(state);

    print("# Starting thermalization... \n");


    #thermalization
    @inbounds for n=1:(nsweeps)*thermfactor
      @inbounds for i=1:(nspins*sweepfactor)
        StateUpdate(state,flips,nflips,1.);
      end;
    end;
    print("# Thermalization done \n");

    #If the number of workers exceeds the number of variational parameters,
    #only a number of workers equals to the number of variational parameters will be used.
    #The exceeding workers will be removed.
    if nprocs() > dim
        println("WARNING: number of workers exceeds number of variational parameters. Only $dim workers will be used.")
        diff = nprocs()-dim
        rmprocs(nprocs()-diff+1:nprocs())
    end
    #Divide the formation of the S-matrix in nprocs() parts, each of them is assigned to a process
    strip = Array{Array{Int64}}(0)
    if nprocs() >= 3
      strip_gen(strip)
    elseif nprocs() <3
        push!(strip,collect(1:dim))
    else error("Number of workers must be a positive integer.")
    end


    #Main loop for the ground state optimization.
    @inbounds for qq=1:n_iterations

        #Sampling of states for nprocs()=1
        if nprocs() == 1
            kk= @parallel (+) for i=1:nprocs()
              sampling(W_RBM,state,nsweeps,i)
          end

           O = Array{Array{Complex{Float64}}}(0)
           for i=1:dim
             push!(O,kk[:,i+1])
           end

           #Formation of the S-matrix
           ss = @parallel (+) for i=1:nprocs()
           S_gen(O,i,strip)
           end
        else
           #Sampling of states for nprocs()>1
           nsweeps_ = Int(round(nsweeps/(nprocs()-1)))     #divide sampling task for each process
           kk= @parallel (+) for i=1:nprocs()-1
             sampling(W_RBM,state,nsweeps_,i)
         end

          O = Array{Array{Complex{Float64}}}(0)
          for i=1:dim
            push!(O,kk[:,i+1])
          end
          #Formation of the S-matrix
          ss = @parallel (+) for i=1:nprocs()-1
          S_gen(O,i,strip)
          end
      end

           ss = ss + eye(dim,dim)*0.0001      #regularization of the S-matrix with a constant shift of the diagonal element of S
           inv_S = inv(ss)                    #inversion of the S-matrix with normal inversion routine (not pseudo-inverse)
           energy = kk[:,1]                   #store of the total energy for each state
           mean_E = mean(energy)              #calculate the mean value of the energy along all the sampled states
           var_E = var(energy)                #Calculate the variance of the energy
           println("iteration $qq- total energy: $(real(mean_E)) +- $(real(var_E)/(sqrt(length(energy))))")

           dF = forces(O,energy)              #calculate the derivative of the energy respect to the variational parameters and assign them to the array dF

    for i=1:dim
      W_RBM[i] = W_RBM[i] - gamma*(inv_S[i,:].'*dF)          #SR optimization rule
    end
    for i in procs()
      fetch(@spawnat i W_to_par(W_RBM,L))                    #update neural network parameters after optimization
    end
  end;

  println("# ground state optimization done");
  println(" ")
end;

#Initialize the variational parameters with small random values
function InitVarPar(nspins,seed)
  srand(seed)
  resize!(a, nspins);
  resize!(b, M);
  resize!(W,nspins);

  for i=1:nhv[2]
      W[i]=typer*zeros(M)
  end
  for i=1:nhv[2]
      a[i] = 0.*typer
  end
  resize!(W_RBM,0)

    for i=1:dim
      push!(W_RBM,rand(-0.01:0.001:0.01) .+ im.*rand(-0.01:0.001:0.01))
    end

  W_to_par(W_RBM,L)

end;

#computes the logarithm of Psi(state')/Psi(state) where
#state' is a state with a certain number of flipped spins
# the vector "flips" contains the sites to be flipped
#look-up tables are used to speed-up the calculation
function Log_wf(state,flips)
      if length(flips) == 0
          return typer*0.;
      end;

      logwf = typer*0.0;

      #Change due to the visible bias
        for i=1:length(flips)
          logwf -= a[flips[i]]*2*state[flips[i]];
        end;

        #Change due to the interaction weights
        for  h=1:M
           thetah=Lt[h];
           thetahp=thetah;

          for i=1:length(flips)
            thetahp -= 2.*state[flips[i]]*W[flips[i]][h];
          end;
          logwf += log(cosh(thetahp))-log(cosh(thetah));
        end;

        return logwf;
end;

function ExpLog(state,flips)
    return exp(Log_wf(state,flips));
end;

#initialization of the look-up tables
function InitLt(state)

  resize!(Lt,M);

  for h=1:M
    Lt[h]=b[h];
    for v=1:nhv[2]
      Lt[h] += state[v]*W[v][h]
    end;
  end;
end;

#updates the look-up tables after spin flips
#the vector "flips" contains the indices of sites to be flipped
function UpdateLt(state,flips)
  if length(flips)==0
    return;
  end;

  for h=1:M
    for i=1:length(flips)
      Lt[h] -= 2*state[flips[i]]*W[flips[i]][h];
    end;
  end;
end;

#total number of spins equal to the number of visible units
function Nspins()
      return nhv[2];
end;

#Generate a spin-flip in the state
function GenSpin(state,flips,nflips,mag0)
     resize!(flips,nflips);

      flips[1]= rand(1:nspins)
      if nflips==2
        flips[2]=rand(1:nspins);
        if !mag0
          return flips[2]!= flips[1];
        else
          return state[flips[2]]!= state[flips[1]];
        end;
      end;

      return true;
end;

#Initializes a random spin state
#if mag0=true, the initial state is prepared with zero total magnetization
function GenRandomState(mag0)
      resize!(state,nspins);
      for i=1:nspins
        state[i]=2*rand(0:1)-1
      end;

    if mag0
        magt=1;
        if nspins%2 == 1
          print("# Error : Cannot initializate a random state with zero magnetization for odd number of spins")
          error();
        end;
        while magt!=0
          magt=0;
          for i=1:nspins
            magt += state[i];
          end;
          if magt>0
            rs = rand(1:nspins);
            while state[rs]<0
              rs= rand(1:nspins);
            end;
            state[rs]=-1;
            magt-=1;

          elseif magt<0
             rs = rand(1:nspins);
            while state[rs]>0
              rs = rand(1:nspins);
            end;
            state[rs]=1;
            magt+=1;

        end;
      end;
    end;

end;


#Generate a new state from the input "state" through a flip of the spins according to the acceptance probability and then update look-up table
function StateUpdate(state::Array{Int64},flips,nflips::Int64,beta::Float64)

    #Picking "nflips" random spins to be flipped
    if GenSpin(state,flips,nflips,mag0)

      #Computing acceptance probability
    acceptance = (norm(ExpLog(state,flips))^2)^beta;

      #Metropolis-Hastings test
      if acceptance>rand()

        #Updating look-up tables in the wave-function
        UpdateLt(state,flips);

        #Moving to the new configuration
        for i=1:length(flips)
          state[flips[i]]*=-1;
        end;

        #AccMov[1]+=1;
      end;
    end;

    #AccMov[2]+=1;
end;

#Measuring the value of the local energy on the current state
function Energy(state,energy)

     en=0.*typer;

    #Finds the non-zero matrix elements of the hamiltonian
    #on the given state
    #i.e. all the state' such that <state'|H|state> = conn(state') \neq 0
    #state' is encoded as the sequence of spin flips to be performed on state
    Pushconn(state,flipsh,conn);

    for i=1:length(flipsh)
      en += ExpLog(state,flipsh[i])*conn[i];
    end;

    push!(energy,en);

end;

#Useful function for the formation of the S-matrix for parallel simulations
function strip_gen(strip::Array{Array{Int64}})
    if nprocs() == 3
        dd = Int(round(dim/3))
        push!(strip,collect(1:dd))
        push!(strip,collect(dd+1:dim))
    end

    dd = Int(round(dim/(2*(nprocs()-1))))

    for i=1:2
      push!(strip,collect(1+(i-1)*dd:i*dd))
    end
    for i=3:nprocs()-2
      push!(strip,collect(1+2*dd+(i-3)*2*dd:(i-2)*2*dd+2*dd))
    end
        push!(strip,collect(1+strip[end][end]:dim))
end

#Form the S-matrix
function S_gen(kk::Array{Array{Complex{Float64}}},y::Int64,strip::Array{Array{Int64}})

       S = typer*zeros(dim,dim)
       RO= real(kk)
       IO= imag(kk)
       meanRO = Array{Complex{Float64}}(dim)
       meanIO = Array{Complex{Float64}}(dim)
       for i=1:dim
           meanRO[i] = mean(RO[i])
           meanIO[i] = mean(IO[i])
       end
       for i in strip[y]
           for j=i:dim
               fR = mean(RO[i].*RO[j]+IO[i].*IO[j])-meanRO[i]*meanRO[j]-meanIO[i]*meanIO[j]
               fI = mean(RO[i].*IO[j]-IO[i].*RO[j])-meanRO[i]*meanIO[j]+meanIO[i]*meanRO[j]
               S[i,j] = fR + im*fI
               S[j,i] = fR -im*fI
             end;
             end;
     RO=0.
     IO=0.
     return S
     end


#Calculate derivatives of the energy respect to the variational parameters
    function forces(O,energy)
        dF= Array{Complex{Float64}}(0)
        mean_E = mean(energy)
        for i=1:dim
            f = mean(energy.*conj(O[i])) - mean(conj(O[i]))*mean_E
            push!(dF,f)
        end
        return dF
    end


#Calculate the logarithm of the wavefunction projected onto the input state
function RBM_wf(state)

  rbm=0.0;
  thetah= 0.0;
  for v=1:nhv[2]
    rbm += a[v]*state[v];
  end;

      for h=1:nhv[1]*nhv[2]
        thetah = b[h];
      for v=1:nhv[2]
          thetah += state[v]*W[v][h];
        end;

        rbm +=log(2cosh(thetah));
      end;

      return exp(rbm);
end;

############## FUNCTIONS FOR DYNAMICS #######################################################################

#Start real-time evolution from initial variational parameters init_par and with "method" integration scheme
function run_dynamics(method,init_par)
    time=collect(0.:step_size:length_int);    #initialized the time-interval of the simulation. Don't touch it!
    println("# Starting dynamics for $n_spins spins and α=$α, $Sym.")
    println("# Time evolution = [0.,$length_int], time-step = $step_size. ")
    println("# Number of sweeps = $nsweeps_d")
    println("# Number of workers = $(nprocs())")
    TT = Array{Complex{Float64}}(dim)
    for i=1:dim
        TT[i] = init_par[i]+0*im
    end
    @time method(rhs_gen,TT,nsweeps_d,time,step_size)
end

#Generate the right hand side of the TDVP equations  d/dt W(t) = -i*S^-1 F(W(t))
function rhs_gen(W_t::Array{Complex{Float64}},t::Float64,nsweeps::Int64,time_flag::Bool)


    nflips=SpinFlips()

    #If the number of workers exceeds the number of variational parameters,
    #only a number of workers equals to the number of variational parameters will be used.
    #The exceeding workers will be removed.
    if nprocs() > dim
        println("WARNING: number of workers exceeds number of viariational parameters. Only $dim workers will be used.")
        diff = nprocs()-dim
        rmprocs(nprocs()-diff+1:nprocs())
    end

    strip = Array{Array{Int64}}(0)
    if nprocs() >= 3
      strip_gen(strip)
  elseif nprocs() <3
      push!(strip,collect(1:dim))
  else error("Number of workers must be a positive integer.")
    end

    #perturbation
    if t >= 0 && t <= 0.2
    @everywhere    pert[1] = true
     else @everywhere pert[1] = false
    end

    for i in procs()
    fetch(@spawnat i W_to_par(W_t,L))
    end

    #initialize a random state
    GenRandomState(mag0);
    resize!(flips,nflips);

    #initializing look-up tables in the wave-function
    InitLt(state);

    #thermalization
    for n=1:nsweeps
      for i=1:(nspins)
        StateUpdate(state,flips,nflips,1.);
      end;
    end;



    #sampling of states if nprocs=1
    if nprocs() == 1
    kk= @parallel (+) for i=1:nprocs()
          sampling(W_t,state,nsweeps,i)
      end

      O = Array{Array{Complex{Float64}}}(0)
      for i=1:dim
        push!(O,kk[:,i+1])
      end

      #Formation of the S-matrix
      s_matrix = @parallel (+) for i=1:nprocs()
        S_gen(O,i,strip)
      end

      #sampling of states if nprocs>1
        else
            nsweeps_ = Int(round(nsweeps/(nprocs()-1)))     #divide sampling task for each process
            kk= @parallel (+) for i=1:nprocs()-1
            sampling(W_t,state,nsweeps_,i)
        end

        O = Array{Array{Complex{Float64}}}(0)
        for i=1:dim
            push!(O,kk[:,i+1])
        end

        #formation of the S-matrix if nprocs()>1
        s_matrix = @parallel (+) for i=1:nprocs()-1
            S_gen(O,i,strip)
        end

    end

      energy = kk[:,1]                      #store the energy for each sampled state
      mean_E = mean(energy)                 #Calculate the mean energy
      var_E = var(energy)                   #Calculate the variance of the energy
      if time_flag == true
      println("time $t- total energy: $(real(mean_E)) +- $(real(var_E)/length(energy))")         #print the energy only in the first step of the Heun scheme
      end
      dF = forces(O,energy)                 #Calculate the derivative of the energy respect to the variational parameters
      dot_W = minres(s_matrix,-im*dF)       #Calculate the rhs of the TDVP equation using minres


      kk=0.
      s_matrix=0.
      O=0.
    return dot_W
end


#Heun integration scheme
function heun(f,W_s::Array{Complex{Float64}},nsweeps::Int64,tspan::Array{Float64},h::Float64)
    push!(W_RBM_t,W_s)

    W_ = Array{Complex{Float64}}(dim)
    W_ .= W_s

    for t in tspan

            l=f(W_RBM_t[end],t,nsweeps,true)
            W_[:] +=  h*(l+ f(W_RBM_t[end]+ h*l,t+h,nsweeps,false))/2
        push!(W_RBM_t, W_[:] )
    end
end
