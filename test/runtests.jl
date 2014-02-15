using Base.Test
require("SDE")
include("../src/misc.jl")
include("../src/quad.jl")

test_type = length(ARGS) == 1 ? ARGS[1] : "TEST"

include("testdiff.jl")
include("testsch.jl")
include("testnonpar.jl")
include("testlinproc.jl")

if test_type == "ALL"
	include("testtc.jl")
end

