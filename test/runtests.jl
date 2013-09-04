using Base.Test
require("SDE")
include("../src/misc.jl")
include("../src/quad.jl")

test_type = length(ARGS) == 1 ? ARGS[1] : "TEST"

include("testdiff.jl")
include("testsch.jl")
include("testlyap.jl")
include("testnonpar.jl")

if test_type == "ALL"
	include("testlinproc.jl")
	include("extendedtestclark.jl")
end

