#define module WvDerivativeTi containing relevant functions
#to calculate the derivative of the RBM wavefunction
#for a given state and for translation invariance symmetry
module WvDerivativeTi


#Return a matrix that is used to sort the N independent weights in the matrix W_ij
function cod(y,L_x,L_y)

  #define vectors for useful storing of indexes
  A = Array{Array{Int64}}(undef,L_x)
  B = Array{Array{Int64}}(undef,L_y)
  C = Array{Array{Int64}}(undef,L_x*L_y)

  for i=1:L_y
    B[i] = collect(1+y*L_x*L_y+(i-1)*L_x:i*L_x+y*L_x*L_y)
  end

  for i=1:L_x
      A[i] = Int64[]
  end

  q = 1
  while q <= L_x
      for i=1:L_y
          append!(A[q],circshift(B[i],-(q-1)))
      end

      q += 1

  end

  for i=1:L_x*L_y

      C[i] = Int64[]
  end

  k = 1
  while k < L_x*L_y
      for h=0:L_y-1
          for s=1:L_x
              append!(C[k],circshift(A[s],-L_x*h))
              k = k + 1
          end
      end
  end

  Z =   Array{Int64}(undef,L_x*L_y,L_x*L_y)

  for i=1:L_x*L_y
      for j=1:L_x*L_y
          Z[i,j] = C[i][j]
      end
  end

  return Z

end


#Define useful objects for the implementation of symmetry
#and efficient calculation of observables for the translation symmetric RBM.
function List(Cod::Array{Int64},ii::Int64,nspins::Int64)

    K = Array{Int64}(undef,nspins)

    for t=1:nspins
        K[t]= findfirst(isequal(ii),Cod[:,t])
    end

    return K

end

#Return a matrix used to efficiently calculate the derivatives
#of the wavefunction for the translation symmetric RBM
function transl_inv_tool(nspins::Int64, alpha::Int64,L_x::Int64,L_y::Int64)

  #useful matrices for filling in the filter matrix
  Cod = Array{Array{Int64}}(undef,alpha)
  ss = Array{Array{Int64}}(undef,alpha)

  #Matrix storing the filter index matrix for transaltion invariance
  filter_ind = Array{Int64}(undef,nspins,nspins*alpha)

  for h=1:alpha
  Cod[h] = cod(h-1,L_x,L_y)
  ss[h] = Cod[h][:,1]
  end

  for h=1:alpha
      for i=1:nspins
          filter_ind[:,i+(h-1)*nspins] = List(Cod[h],ss[h][i],nspins)
      end
  end

  return filter_ind

end

#Return a vector used to efficiently calculate the derivatives
#of the wavefunction for the translation symmetric RBM
function ti_tool(nspins::Int64, L_x::Int64,L_y::Int64,ind_h::Int64,ind_i::Int64)

  Cod = cod(ind_h - 1,L_x,L_y)
  ss = Cod[:,1]

  return List(Cod,ss[ind_i],nspins)

end

#Calculate the derivative of the wavefunction respect to the variational parameters for a given state
function Var_derivatives(state::Array{Int64,1},ENWfDer::Array{Complex{Float64},2},Lt::Array{Complex{Float64},1},FilterInd::Array{Int64,2},nspins::Int64,alpha::Int64, ind::Int64)

    #Fill ENWfDer with the derivative of the wavefunction w.r.t. b and W
    k = 1
    for h=1:alpha
        r_ = (h-1)*nspins; @views t_ = tanh.(Lt[(r_+1):(h*nspins)])
        @views    ENWfDer[ind,1 + h] = sum(t_)
        for i=1:nspins
                @views    lt = state[FilterInd[:,i+r_]]'*t_
                    ENWfDer[ind,1 + alpha + k] = lt
                    k += 1
        end
    end

end

end
