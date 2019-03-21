const dim =M*nhv[2] + M + nhv[2]        #number of variational parameters

#calculate the derivative of the wavefunction respect to the variational parameters for a given state
function Var_derivatives(state,O_S,O_b,O_w)

  for i=1:nhv[2]
  push!(O_S[i],state[i])
end
  for h=1:M
  push!(O_b[h],tanh(Lt[h]))
end

for h=1:M
  for i=1:nhv[2]
    push!(O_w[i,h],state[i]*tanh(Lt[h]))
  end
end
end

#Sample the states and calculate the energy and the derivative of the wavefunction for all the states
function sampling(W_RBM,st::Array{Int64},nsweeps::Int64,y::Int64;thermfactor=0.1,sweepfactor=1)
    srand(2245*y)
    nflips=SpinFlips()
    flips_p =  Array{Int64}(2)
    state_p = st

    energy = Array{Complex{Float64}}(0)
    O_S = Array{Array{Int64}}(nhv[2])
    O_b = Array{Array{Complex{Float64}}}(M)
    O_w = Array{Array{Complex{Float64}}}(M*nhv[2])
    O_w = reshape(O_w,nhv[2],M)

   resize!(flips_p,nflips);

   #initializing look-up tables in the wave-function
   InitLt(state_p);

   ReStateUpdate(O_S,O_b,O_w)

   #sequence of sweeps

   @inbounds for n=1:nsweeps
       @inbounds    for i=1:(nspins*sweepfactor)
           StateUpdate(state_p,flips_p,nflips,1.);
       end;

       Energy(state_p,energy)

       Var_derivatives(state_p,O_S,O_b,O_w)

   end;
  state[:]=state_p

  if nprocs() == 1
        q=typer*zeros(nprocs()*nsweeps,1+dim)
  else  q=typer*zeros((nprocs()-1)*nsweeps,1+dim)
  end

 reshape(O_w,M*nhv[2])
 for i=1:nsweeps
     q[nsweeps*(y-1)+i,1]=energy[i]
 end
 for i=1:nsweeps
     for d=1:nhv[2]
         q[nsweeps*(y-1)+i,d+1] = O_S[d][i]
     end
 for j=1:M
     q[nsweeps*(y-1)+i,j+1+nhv[2]] = O_b[j][i]
 end
 for k=1:M*nhv[2]
     q[nsweeps*(y-1)+i,k+1+M+nhv[2]]=O_w[k][i]
 end
 end

 return q
end

#update the neural network parameters
function W_to_par(W_RBM,L)
  a[:] = W_RBM[1:nhv[2]]
  b[:] = W_RBM[(nhv[2]+1):(nhv[2]+M)]
  q=1
  while (q+M+nhv[2]) <= dim
    for j=1:M
      for i=1:nhv[2]
        W[i][j] = W_RBM[M+nhv[2]+q]
        q +=1
      end;
    end;

end;
end;

#useful function for gc()
function ReStateUpdate(O_S::Array{Array{Int64}},O_b::Array{Array{Complex{Float64}}},O_w::Array{Array{Complex{Float64}}})
  for i=1:nhv[2]
    O_S[i]=Int64[]
  end;
  for h=1:M
    O_b[h]=Complex{Float64}[]
    for i=1:nhv[2]
    O_w[i,h]=Complex{Float64}[]
  end
  end

end;
