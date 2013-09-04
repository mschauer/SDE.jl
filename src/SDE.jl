#include(Pkg.dir("SDE")*"/src/Diffusion.jl")
#include(Pkg.dir("SDE")*"/src/Schauder.jl")
#include(Pkg.dir("SDE")*"/src/NonparBayes.jl")
#include(Pkg.dir("SDE")*"/src/Randm.jl") 
#include(Pkg.dir("SDE")*"/src/LinProc.jl")
#include(Pkg.dir("SDE")*"/src/Lyap.jl")

#include relative to location of SDE.jl, not to current pwd()
include("Diffusion.jl")
include("Schauder.jl")
include("NonparBayes.jl")
include("Randm.jl") 
include("LinProc.jl")
include("Lyap.jl")


#placeholder
module SDE

end
