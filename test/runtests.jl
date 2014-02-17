using Base.Test
require("SDE")
#is called with arg "travis" from travis
test_type = length(ARGS) == 1 ? ARGS[1] : ""

include("testsde.jl")
include("testlinproc.jl")

if test_type == "all"
	include("testtc.jl")
end

include("testsch.jl")
include("testnpbayes.jl")
