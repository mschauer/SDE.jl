using Base.Test
require("SDE")
#include("../src/misc.jl")
#include("../src/quad.jl")

test_type = length(ARGS) == 1 ? ARGS[1] : ""

include("testsde.jl")
include("testsch.jl")
include("testlinproc.jl")

if test_type == "all"
	include("testtc.jl")
end

include("testnonpar.jl")
