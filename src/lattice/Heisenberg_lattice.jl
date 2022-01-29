#THis module defines the lattice geometry, namely the bonds
#in the lattice with or without periodic boundary condition (Pbc)
#for the Heisenberg model in a L_x x L_y lattice

module Lattice

#Initialization of the lattice (store all the bonds in the Heisenberg Hamiltonian)
function InitLattice(L_x::Int64, L_y::Int64, nspins::Int64)

    #Define arrays that will be filled with nearest neigbourng bonds
    v_bonds= Array{Array{Int64}}(undef,0)
    h_bonds= Array{Array{Int64}}(undef,0)

    #Set boundary conditions
    if L_x > 2 && L_y > 2
        pbc_v = true; pbc_h = true
    elseif L_x <= 2 && L_y > 2
        pbc_v = true; pbc_h = false
    elseif L_x > 2 && L_y == 2
        pbc_v = false; pbc_h = true
    else
        pbc_v = false; pbc_h = false
    end


    #Check if the number of spins is compatible with the given lattice size
    if (L_x*L_y != nspins)
        println("Error , the number of spins is not compabitle with lattice length ")
        error();
    end;

    nn = Array{Array}(undef,nspins)
    for i=1:nspins
        nn[i]= Array{Int64}(undef,4)
    end;

    #Defining the nearest-neighbors
    for i=1:nspins

        if pbc_h
            nn[i][1] = PbcH(L_x,i-1,i)
        else nn[i][1] = PbcH(L_x,i-1,i)
        end

        if pbc_h
            nn[i][2] = PbcH(L_x,i+1,i)
        else nn[i][2] = PbcH(L_x,i+1,i)
       end

       if pbc_v
           nn[i][3] = PbcV(nspins,i-L_x)
       else nn[i][3] = V(nspins,i-L_x)
       end

       if pbc_v
           nn[i][4] = PbcV(nspins,i+L_x)
       else nn[i][4]= V(nspins,i+L_x)
       end
   end

 #Push nearest neighbourong bonds along y into the array v_bonds
   for  i=1:nspins
       for  k=3:4
           j = nn[i][k]
           if i<j
               push!(v_bonds,[i;j])
           end
       end
   end

  #Push nearest neighbourong bonds along x into the array h_bonds
   for  i=1:nspins
       for  k=1:2
           j = nn[i][k];
           if i<j
               push!(h_bonds,[i;j])
           end
       end
   end

   #We only keep inequivalent bonds
    bondsH = unique(h_bonds)
    bondsV = unique(v_bonds)

   return bondsH, bondsV

end


#Small functions to set up the lattice Horizontal (x) Pbc
function PbcH(L_x,nn,s)
    if (s-1)%L_x==0 && nn== (s-1)
      return (s+L_x-1);
  elseif (s)%L_x==0 && nn==(s+1)
      return (s-L_x+1);
    else
      return nn;
  end;
end;


#Small functions to set up the lattice Vertical (y) Pbc
function PbcV(nspins,nn)
    if nn>nspins
      return (nn-nspins);
  elseif nn<=0
      return (nspins+nn);
    else
      return nn;
  end;
end;


#Small functions to set up the lattice Horizontal (x) without Pbc
function H(L_x,nn,s)
    if s%L_x==0 && nn==(s-1)
      return -1
  elseif (s+1)%L_x==0 && nn==(s+1)
      return -1
    else
      return nn
  end;
end;


#Small functions to set up the lattice Vertical (y) without Pbc
function V(nspins,nn)
    if nn>nspins
      return -1;
    elseif nn<0
      return -1;
    else
      return nn;
  end;
end;

end
