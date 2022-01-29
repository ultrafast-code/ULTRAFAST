#Define function to calculate gradient of energy
function Egrad(dim::Int64)

    dE = Array{Complex{Float64}}(undef,0)

    for i=1:dim
    @views     push!(dE,cov(ENWfDer[:,1],ENWfDer[:,i+1],corrected=false))
    end

    return dE

end
