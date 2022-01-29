# MyULTRAFAST


MyULTRAFAST is an open-source code for the simulation of ground state and dynamic
properties of quantum spin systems.

It is derived from the open-source code ULTRAFAST available at
        
        https://github.com/ultrafast-code/ULTRAFAST 

written by G. Fabiani and J.H. Mentink for research purpose. Any use of MyULTRAFAST for scientific research has to provide citation to ULTRAFAST.

MyULTRAFAST exploits the learning power of artificial neural networks, in particular of the restricted Boltzmann machine architecture, to simulate quantum states of spin-1/2 Hamiltonians in one and two dimensions. By using the high-level parallelization routines provided by Julia, MyULTRAFAST offers the possibility to fully parallelize the simulation over the available cores, allowing to simulate lattice sizes as large as 30x30 in modest CPU time. It also offers full customization of the model Hamiltonian, the possibility to adapt the code to any lattice geometry and implement any dynamic term.


To run MyULTRAFAST, first install Julia, version 1.1.1 or higher. Then install the code by 
downloading the tar file given in the Supplementary Material and untarring it in 
a suitable working directory. Julia can be run either from an interactive session
Read-Eval-Print Loop ("REPL")  or from the command line. To execute ULTRAFAST
in the REPL, double-click the Julia executable and type

julia> include("../src/main/main.jl")

From the command line, open a terminal and type

$ path/to/julia   path/to/src/main/main.jl

To get acquainted with MyULTRAFAST, we refer to the tutorials in the Example folder.
