const dim = nhv[1] + M                                #dimension of the variational manifold
const lista = Array{Int64}(nhv[2],nhv[1]*nhv[2])      #Useful quantity in the optimization
const Cod = Array{Array{Int64}}(nhv[1])               #Useful quantity in the optimization
const ss = Array{Array{Int64}}(nhv[1])                #Useful quantity in the optimization

#Define useful objects for the implementation of symmetry and efficient calculation of observables for the translation symmetric wavefunction.
function transl_inv_tool()
  for h=1:nhv[1]
  Cod[h] = cod(h-1,L)
  ss[h]=Cod[h][:,1]
  end

  for h=1:nhv[1]
  for i=1:nhv[2]
    lista[:,i+(h-1)*nhv[2]] = List(Cod[h],ss[h][i],L)
  end
  end
end

#Calculate the derivative of the wavefunction respect to the variational parameters for a given state
function Var_derivatives(lis::Array{Int64},state,O_b,O_w)
    @inbounds    for h=1:nhv[1]
                    push!(O_b[h],sum(tanh.(Lt[((h-1)*nhv[2]+1):(h*nhv[2])])))
                 end

       for h=1:nhv[1]
        for i=1:nhv[2]
                       lt = state[lis[:,i+(h-1)*nhv[2]]]'*tanh.(Lt[1+(h-1)*nhv[2]:h*nhv[2]])
                       push!(O_w[i+(h-1)*nhv[2]],lt)
        end
       end
end

#Sample nsweeps states and calculate energy and wavefunction derivatives for all the states
function sampling(W_RBM,st::Array{Int64},nsweeps::Int64,y::Int64;thermfactor=0.1,sweepfactor=1)
       srand(100*y)
       nflips=SpinFlips()
       flips_p =  Array{Int64}(40)
       state_p = st

       energy = Array{Complex{Float64}}(0)
       O_b = Array{Array{Complex{Float64}}}(nhv[1])
       O_w = Array{Array{Complex{Float64}}}(M)

       resize!(flips_p,nflips);

       #initializing look-up tables in the wave-function
       InitLt(state_p);

       ReStateUpdate(O_b,O_w)

       #sequence of sweeps

       @inbounds for n=1:nsweeps
           @inbounds    for i=1:(nspins*sweepfactor)
               StateUpdate(state_p,flips_p,nflips,1.);
           end;

           Energy(state_p,energy)

           Var_derivatives(lista,state_p,O_b,O_w)

       end;

     if nprocs() == 1
           q=typer*zeros(nprocs()*nsweeps,1+dim)
     else  q=typer*zeros((nprocs()-1)*nsweeps,1+dim)
     end

     for i=1:nsweeps
         q[nsweeps*(y-1)+i,1]=energy[i]
     end
     for i=1:nsweeps
      for j=1:nhv[1]
         q[nsweeps*(y-1)+i,j+1] = O_b[j][i]
     end
     for k=1:M
         q[nsweeps*(y-1)+i,k+nhv[1]+1] = O_w[k][i]
     end
     end

     return q
end


#Return a matrix that is used in W_to_par() to sort the N independent weights in the matrix W_ij
function cod(y,L)
  B = Array{Array{Int64}}(L)
  A = Array{Array{Int64}}(L)
  string = Array{Array{Int64}}(L*L)
  for i=1:L
    B[i]=collect(1+y*L^2+(i-1)*L:i*L+y*L^2)
    A[i]=Int64[]
  end
  q=1
  while q<=L
  for i=1:L
    append!(A[q],circshift(B[i],-(q-1)))
  end

  q +=1
end

for i=1:L^2
string[i]=Int64[]
end
  k=1
  while k < L^2
  for s=1:L
  for h=0:L-1
    append!(string[k],circshift(A[s],-L*h))
    k = k+ 1
  end

end
end
Z=Array{Int64}(L^2,L^2)
for i=1:L^2
  for j=1:L^2
      Z[i,j]= string[i][j]
  end
end
return Z
end

#Return a matrix used to efficiently calculate the derivatives of the wavefunction for the translation symmetric wavefunction
function List(Cod::Array{Int64},ii::Int64,L::Int64)
    K = Array{Int64}(L^2)
    for t=1:L^2
        K[t]= findfirst(Cod[:,t],ii)
    end
return K
end

#Update the neural network parameters after the SR update of the variational parameters, assigning the variational parameters to b_j and W_ij
function W_to_par(W_t,L::Int64)
  C = Array{Complex{Float64}}(0)
  W_pre = Array{Array{Int64}}(nhv[1])

  for h=1:nhv[1]
    W_pre[h] = cod(h-1,L)
  end

  for i=1:nhv[1]
    append!(C,W_t[i].*ones(nhv[2]))
  end

  W_s=W_t[nhv[1]+1:end]
  b[:] = C

  for i=1:L^2
    W[i]= Complex{Float64}[]
    for h=1:nhv[1]
      append!(W[i],W_s[W_pre[h][:,i]])
    end
  end
end

#gc function for the derivative of the wavefunction O_b, O_w
function ReStateUpdate(O_b::Array{Array{Complex{Float64}}},O_w::Array{Array{Complex{Float64}}})
  for h=1:nhv[1]
    O_b[h]=Complex{Float64}[]
  end
  for i=1:nhv[1]*nhv[2]
    O_w[i]=Complex{Float64}[]
  end
end;


transl_inv_tool()
