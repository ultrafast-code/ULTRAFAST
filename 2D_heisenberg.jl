
const bonds= Array{Array{Int64}}(0)
const pert_bonds= Array{Array{Int64}}(0)
const flipsh= Array{Array{Int64}}(0)
const conn = Array{Complex{Float64}}(1)

function Heisenberg2d(nspins,coupling,pbc)
  InitLattice()
end;

#Initialization of the lattice (store all the bonds in the Hamiltonian)
function InitLattice()
 if (L*L != nspins)
   println("Error , the number of spins is not compabitle with a square lattice ")
   error();
 end;

 nn=Array{Array}(nspins)
 for i=1:nspins
 nn[i]= Array{Int64}(4)
 end;

 #Defining the nearest-neighbors
 for i=1:(nspins)
   if pbc
    nn[i][1] = pbcx_bonds(i-1,i)
  else nn[i][1]= H(i-1,i)
     end;
   if pbc
      nn[i][2] = pbcx_bonds(i+1,i)
    else nn[i][2]= H(i+1,i)
       end;
    if pbc
      nn[i][3] = pbcy_bonds(i-L)
  else nn[i][3]= V(i-L)
         end;
    if pbc
        nn[i][4] = pbcy_bonds(i+L)
    else nn[i][4]= V(i+L)
           end;
 end;

 for  i=1:nspins
   for  k=1:4
      j=nn[i][k];
     if i<j
       push!(bonds,[i;j])
     end;
   end;
 end;
 for  i=1:nspins
   for  k=3:4
      j=nn[i][k];
     if i<j
       push!(pert_bonds,[i;j])
     end;
   end;
 end;

end;



#Finds the non-zero matrix elements of the hamiltonian on the given state
#i.e. all the state' such that <state'|H|state> = conn(state') \neq 0
#state' is encoded as the sequence of spin flips to be performed on state
function Pushconn(state,flipsh,conn)

  resize!(conn,1);
  resize!(flipsh,1);
  flipsh[1]= Int64[];
  #computing interaction part Sz*Sz
  conn[1]=0.*typer;

    if nspins == 4
        count = 2
    else count = 1
    end

    #Contribution to the energy of the zz component of the Heisenberg Hamiltonian
    for i=1:count:length(bonds)
      conn[1] += state[bonds[i][1]]*state[bonds[i][2]];
    end;
    conn[1]*= jz;


    #Looks for possible spin flips
    for i=1:count:length(bonds)
      const  si=bonds[i][1];
      const  sj=bonds[i][2];

      #Contribution of the off-diagonal terms to the total energy.

      if state[si]!=state[sj]
        push!(flipsh,[si;sj]);
        push!(conn,-2.0);
      end;
    end;

    #Contribution to the total energy of the perturbation term due to the perturbation of the exchange interaction along vertical bonds
    if pert[1]
        for i=1:1:length(pert_bonds)
          conn[1] += coupling[1]*state[pert_bonds[i][1]]*state[pert_bonds[i][2]];
        end;

        #Looks for possible spin flips
        for i=1:count:length(pert_bonds)
          const  si=pert_bonds[i][1];
          const  sj=pert_bonds[i][2];

          if state[si]!=state[sj]
            push!(conn,-2.0*coupling[1]);
            push!(flipsh,[si;sj]);
          end;
        end;
    end
end
#Select the number of spin flips
function SpinFlips()
    return n_flips;
  end;


#Small functions to set up the lattice Horizontal Pbc
  function pbcx_bonds(nn,s)
    if (s-1)%L==0 && nn==(s-1)
      return (s+L-1);
  elseif (s)%L==0 && nn==(s+1)
      return (s-L+1);
    else
      return nn;
  end;
end;
  #Vertical Pbc
  function pbcy_bonds(nn)
    if nn>nspins
      return (nn-nspins);
  elseif nn<=0
      return (nspins+nn);
    else
      return nn;
  end;
end;
  #Horizontal without pbc
   function H(nn,s)
    if s%L==0 && nn==(s-1)
      return -1;
  elseif (s+1)%L==0 && nn==(s+1)
      return -1;
    else
      return nn-1;
  end;
end;

  #Vertical without pbc
  function V(nn)
    if nn>=nspins
      return -1;
    elseif nn<0
      return -1;
    else
      return nn;
  end;
end;


Heisenberg2d(nspins,coupling,pbc)
