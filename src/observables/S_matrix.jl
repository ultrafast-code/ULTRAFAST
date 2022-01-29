#Define functions for the calculation and formation of the covariance matrix S

#Compute the element of the S-matrix and store them into a vector
function S_gen(dim::Int64)

       S_vec = Complex{Float64}[]

       for i=1:dim
           for j=i:dim
      @views   fr = cov(ENWfDer[:,j+1],ENWfDer[:,i+1],corrected=false)
               push!(S_vec,fr)
             end
       end

     return S_vec

end

#Form and return the S-matrix
function Gen_S(S_vec::Array{Complex{Float64},1})

    S_mat = zeros(Complex{Float64},dim,dim)

    k = 1
    for i=1:dim
        for j=i:dim
            @views S_mat[i,j] = S_vec[k]
            @views S_mat[j,i] = conj(S_vec[k])
            k +=1
        end
    end

    return S_mat

end
