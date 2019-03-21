#Set Neural Network
 n_spins = 16               #number of spins (visible units)
 α =  4                     #number hidden units/visible units
 const pbc = true           #if true ==> periodic boundary conndition

#Set symmetries
#Uncomment the symmetry you want to employ
#Sym = "No symmetry"
Sym = "Translation symmetry"

mag0 = true                 #if true ==> zero-magnetization sector
const n_flips = 2           #spin flips in the monte carlo sampling. Set n flips=2 for mag0=true
