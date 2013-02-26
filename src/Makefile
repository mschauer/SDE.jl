all: ../doc/diffusion.rst ../doc/schauder.rst ../doc/sdenonparbayes.rst

../doc/diffusion.rst: diffusion.jl
	sed -ne "s/^#% \? \?//p" diffusion.jl > ../doc/diffusion.rst
	
../doc/schauder.rst: schauder.jl
	sed -ne "s/^#% \? \?//p" schauder.jl > ../doc/schauder.rst
	
../doc/sdenonparbayes.rst: sdenonparbayes.jl
	sed -ne "s/^#% \? \?//p" sdenonparbayes.jl > ../doc/sdenonparbayes.rst
