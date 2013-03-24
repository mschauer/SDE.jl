all: ../doc/diffusion.rst ../doc/schauder.rst ../doc/sdenonparbayes.rst ../doc/randm.rst ../doc/lyap.rst ../doc/linproc.rst

../doc/diffusion.rst: diffusion.jl
	sed -ne "s/^#% \? \?//p" diffusion.jl > ../doc/diffusion.rst
	
../doc/schauder.rst: schauder.jl
	sed -ne "s/^#% \? \?//p" schauder.jl > ../doc/schauder.rst
	
../doc/sdenonparbayes.rst: sdenonparbayes.jl
	sed -ne "s/^#% \? \?//p" sdenonparbayes.jl > ../doc/sdenonparbayes.rst

../doc/randm.rst: randm.jl
	sed -ne "s/^#% \? \?//p" randm.jl > ../doc/randm.rst

../doc/lyap.rst: lyap.jl
	sed -ne "s/^#% \? \?//p" lyap.jl > ../doc/lyap.rst

../doc/linproc.rst: linproc.jl
	sed -ne "s/^#% \? \?//p" linproc.jl > ../doc/linproc.rst

