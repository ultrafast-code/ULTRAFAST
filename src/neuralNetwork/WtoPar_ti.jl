#Define function to convert independent variational parameters
#to network parameters a, b and W with translation symmetry
function WToPar(RBMPar::RBM_par, WRBM::Array{Complex{Float64},1})

  C = Array{Complex{Float64}}(undef,0)
  W_pre = Array{Array{Int64}}(undef,alpha)

  for h=1:alpha
    W_pre[h] = WvDerivativeTi.cod(h-1,L_x,L_y)
  end

  for i=1:alpha
    append!(C,WRBM[i].*ones(nspins))
  end

  W_s = WRBM[alpha+1:end]
  RBMPar.b[:] = C

  for i=1:nspins
    for h=1:alpha
       RBMPar.W[i,1+nspins*(h-1):h*nspins]=W_s[W_pre[h][:,i]]
    end
  end

end

#Define function to convert network parameters a,b and W
#to independent variational parameters with translation symmetry
function ParToW(RBMPar::RBM_par,W_RBM::Array{Complex{Float64},1})

    for h=1:alpha
        W_RBM[h] = RBMPar.b[1 + (h-1)*nspins]
    end

    W_RBM[alpha+1:end] = RBMPar.W[1,:]

end
