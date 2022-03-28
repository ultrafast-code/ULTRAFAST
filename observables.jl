store_latt = Int64[]
const sub_A = Array{Int64}(0)           #store site indices  in sublattice A
const sub_B = Array{Int64}(0)           #store site indices in sublattice B
for i=1:L
    for k=1:Int(L/2)
        push!(store_latt,(-1)^(i)*(i+2*(k-1)*L))
        push!(store_latt,-(-1)^(i)*(i+L+2*(k-1)*L))
    end
end
#deleteat!(store_latt,1)
for i=1:length(store_latt)
    if store_latt[i] < 0
        push!(sub_A,-store_latt[i])
    else push!(sub_B,store_latt[i])
    end
end

#calculate the energy and spin correlation for spin i_ and j_ for a given state
function Measure(energy_,corr_z,corr_x,corr_y,i_::Int64,j_::Int64,gs::Bool)

   s_x = 0.*typer
   s_y = 0.*typer
   en=0.*typer

   if gs == true
   Pushconn(state,flipsh,conn);

   for i=1:length(flipsh)
     en += ExpLog(state,flipsh[i])*conn[i];
   end;
        push!(energy_,en);
   end

   if (i_ in sub_A && j_ in sub_B) || (i_ in sub_B && j_ in sub_A)

       EL_store = ExpLog(state,[i_;j_])
       s_x -= EL_store
       if state[i_] == state[j_]
           s_y +=  EL_store
       else s_y -= EL_store
       end

    push!(corr_x, (s_x))
    push!(corr_y, (s_y))
    push!(corr_z, state[i_]*state[j_])

    else
        EL_store = ExpLog(state,[i_;j_])
        s_x = EL_store
        if state[i_] == state[j_]
            s_y -=  EL_store
        else s_y += EL_store
        end

         push!(corr_x, (s_x))
         push!(corr_y, (s_y))
         push!(corr_z, state[i_]*state[j_])
    end
end;

#sample nsweeps states and apply measure() at each state through Measure_par(). Outputs the value of spin-spin correlation
function spincorr(W_t,nsweeps::Int64,t::Float64,i_::Int64,j_::Int64;thermfactor=0.1,sweepfactor=1)

    nflips=SpinFlips()
    if (t >= 0 && t <= 0.2)
    @everywhere    pert[1] = true
    else @everywhere pert[1] = false
    end

    for i in procs()
    fetch(@spawnat i W_to_par(W_t,L))
    end

    GenRandomState(mag0);

    resize!(flips,nflips);

    #initializing look-up tables in the wave-function
    InitLt(state);

    #thermalization
    for n=1:nsweeps*thermfactor
      for i=1:(nspins*sweepfactor)
       StateUpdate(state,flips,nflips,1.);
      end;
    end;

    if nprocs() == 1
        kk =  Measure_par(W_t,state,nsweeps,1,i_,j_)
    else
        nsweeps_ = Int(round(nsweeps/(nprocs()-1)))     #divide sampling task for each process
        kk= @parallel (+) for i=1:nprocs()-1
            Measure_par(W_t,state,nsweeps_,i,i_,j_)
        end
    end
    return (mean(kk[:,1]) + mean(kk[:,2]) + mean(kk[:,3]))
end

#function that fills the array for the correlations along x,y,z with multiple processes
function Measure_par(W_t,st::Array{Int64},nsweeps::Int64,y::Int64,i_::Int64,j_::Int64;sweepfactor=1)
   srand(2224*y)
   nflips=SpinFlips()
   flips_p =  Array{Int64}(40)
   state_p = st

   corr_x = Array{Complex{Float64}}(0)
   corr_y = Array{Complex{Float64}}(0)
   corr_z = Array{Complex{Float64}}(0)
   energy_ = Array{Complex{Float64}}(0)
   #initializing look-up tables in the wave-function
   resize!(flips_p,nflips);

   #initializing look-up tables in the wave-function
   InitLt(state_p);

   #sequence of sweeps

    for n=1:nsweeps
        for i=1:(nspins*sweepfactor)
            StateUpdate(state_p,flips_p,nflips,1.);
        end;
        Measure(energy_,corr_z,corr_x,corr_y,i_,j_,true)
    end;

    if nprocs() == 1
        q=typer*zeros(nprocs()*nsweeps,3)
    else
      q=typer*zeros((nprocs()-1)*nsweeps,3)
    end

    for i=1:nsweeps
        q[nsweeps*(y-1)+i,1]=corr_z[i]
        q[nsweeps*(y-1)+i,2]=corr_x[i]
        q[nsweeps*(y-1)+i,3]=corr_y[i]
    end

 return q
end

#calculate the time-evolution of spin-spin correlation and return it
function spincorr_d(W_RBM_t::Array{Array{Complex{Float64}}},Nsample::Int64,step_size::Float64,i_,j_)

    if !isa(i_,Int64) || !isa(j_,Int64) || i_ > nhv[2] || j_ > nhv[2] || i_ <= 0 || j_ <= 0
        error("spin index i must be an integer and 1 <= i <= $(nhv[2]) ")
    end
    println("Calculating time-evolution spin correlation...")
    t=Array{Complex{Float64}}(0)

        for i=1:1:length(W_RBM_t)
            s = spincorr(W_RBM_t[i],Nsample,step_size*(i-1),i_,j_)
            push!(t, s)
        end

    return t
end

#generate nsweeps samples and for each of them calculate the energy and spin-spin correlation, which are stored in arrays.
#print_obs() then returns the ground state energy and spin-spin correlation.
function GS_obs(W_RBM,nsweeps::Int64,i_,j_;thermfactor::Float64=0.1,sweepfactor::Int64=1)

    if !isa(i_,Int64) || !isa(j_,Int64) || i_ > nhv[2] || j_ > nhv[2] || i_ <= 0 || j_ <= 0
        error("spin index i must be an integer and 1 <= i <= $(nhv[2]) ")
    end
      println("# Calculating ground state observables...")
      energy_ =  Array{Complex{Float64}}(0)
      corr_x = Array{Complex{Float64}}(0)
      corr_y = Array{Complex{Float64}}(0)
      corr_z = Array{Complex{Float64}}(0)

      W_to_par(W_RBM,L)

  nflips=SpinFlips()

  GenRandomState(mag0);

  resize!(flips,nflips);

  #initializing look-up tables in the wave-function
  InitLt(state);

  #thermalization
 for n=1:(nsweeps*thermfactor)
    for i=1:(nspins*sweepfactor)
      StateUpdate(state,flips,nflips,1.);
    end;
  end;

  resize!(flips,nflips);

  #initializing look-up tables in the wave-function
  InitLt(state);

  #sequence of sweeps

   for n=1:nsweeps
     for i=1:(nspins*sweepfactor)
            StateUpdate(state,flips,nflips,1.);
     end;

        Measure(energy_,corr_z,corr_x,corr_y,i_,j_,true)
  end;
        tot_corr = corr_z + corr_x + corr_y

        println("# Energy per spin:  $(print_obs(energy_,true))")
        println("# <S_m  S_n>:       $(print_obs(tot_corr,false))")
        println(" ")
end;

#calculate the mean energy and spin-spin correlation and the respective errors
function print_obs(obs,EN::Bool)
  nblocks=200;

  blocksize=Int(floor(length(obs)/nblocks));

  mean_obs=0;
  mean_err=0;

   mean_obs_unb=0;
   mean_err_unb=0;

  for i=0:(nblocks-1)
     obsblock=0;
    for j=(i*blocksize):((i+1)*blocksize-1)
      obsblock+=real(obs[j+1]);
      if j>length(obs)
        error();
      end;

      delta = real(obs[j+1])-mean_obs_unb;
      mean_obs_unb += delta/(j+1);
      delta2 =real(obs[j+1])-mean_obs_unb;
      mean_err_unb+=delta*delta2;
    end;
    obsblock /= blocksize;
    delta = obsblock-mean_obs;
    mean_obs += delta/(i+1);
    delta2 = obsblock-mean_obs;
    mean_err += delta*delta2;
  end;

  mean_err /= (nblocks-1);
  mean_err_unb /= (nblocks*blocksize-1);

  estav = mean_obs/nspins;
 esterror = sqrt(mean_err/nblocks)/nspins;

  ndigits = log10(esterror);
  if ndigits<0
    ndigits =- ndigits+2;
  else
    ndigits=0;
  end;

  if EN == true
      return "$(estav/4) +- $(esterror/4)"
  else return "$(estav*nspins/4) +- $(esterror*nspins/4)"

 end;

end
