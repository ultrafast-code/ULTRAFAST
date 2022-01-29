#Define useful functions to handle parallelization

#Return the number of effective workers
#If more than 1 worker is present Julia will use nprocs()-1 workers
function NumberWorkers()
    if nprocs() >= 3
        return nprocs() - 1
    elseif nprocs() <3
        return 1
    else error("Number of workers must be a positive integer.")
    end
end



#Divide the total number of sample in nprocs() equal parts
#If more than 1 worker is present Julia will use nprocs()-1 workers
function partition_sample(nsample::Int64)
    if nprocs() == 1
        return nsample
    else
        return Int(round(nsample/(nprocs()-1)))
    end
end



#Divide number of entris of S-matrix into nprocs()-1 parts
#and sort them into a given array strip
#dim=number of variational parameters
function GenerateStrips(strip::Array{Array{Int64}}, dim::Int64)

    if nprocs() == 3
        dd = Int(round(dim/2))
        push!(strip,collect(1:dd))
        push!(strip,collect(dd+1:dim))
        return
    end

    dd = Int(round(dim/(2*(nprocs()-1))))

    for i=1:2
      push!(strip,collect(1+(i-1)*dd:i*dd))
    end

    for i=3:nprocs()-2
      push!(strip,collect(1+2*dd+(i-3)*2*dd:(i-2)*2*dd + 2*dd))
    end

    push!(strip,collect(1+strip[end][end]:dim))

end
