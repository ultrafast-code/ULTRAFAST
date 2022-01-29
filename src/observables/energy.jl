#define functions to calculate the energy with Monte Carlo sampling

#calculate scalar product S_i \cdot S_j of a bond [i,j]
function sp_bond(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1},bond::Array{Int64,1})

   if state[bond[1]] == state[bond[2]]
        return 1.0
   else return -1.0 - 2.0*ExpLog(RBMPar,Lt,state,bond)
   end

end


#Finds the non-zero matrix elements of the hamiltonian on the given state
#i.e. all the state' such that <state'|H|state> = conn(state') \neq 0
#state' is encoded as the sequence of spin flips to be performed on state
function Pushconn(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1})

    #initialize array of S_i \cdot S_j between neighbouring spin along x
    en_x = Complex{Float64}[]

    #initialize array of S_i \cdot S_j between neighbouring spin along y
    en_y = Complex{Float64}[]

    #fill en_x and en_y with S_i \cdot S_j on all corresponding bonds
    for i=1:length(bonds_x)
    @views    push!(en_x,sp_bond(RBMPar,Lt,state,bonds_x[i]))
    end

    for i=1:length(bonds_y)
    @views    push!(en_y,sp_bond(RBMPar,Lt,state,bonds_y[i]))
    end

    return en_x, en_y
end


#Measuring the value of the local energy on the current state
#and fill it in the given vector energy
function Energy(RBMPar::RBM_par,Lt::Array{Complex{Float64},1},state::Array{Int64,1},energy::Array{Complex{Float64}},J::Array{Float64,1})

    #calculate S_i \cdot S_j for bonds along x and y
    sp_x, sp_y = Pushconn(RBMPar,Lt,state)

    #calculate the energy by multiplying S_i \cdot S_j with corresponding J
    @views en = J[1]*sum(sp_x) + J[2]*sum(sp_y)

    push!(energy,en)

end;
