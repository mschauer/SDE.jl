all: ../doc/diffusion.rst ../doc/schauder.rst ../doc/nonparbayes.rst ../doc/randm.rst ../doc/lyap.rst ../doc/linproc.rst

../doc/diffusion.rst: Diffusion.jl
	sed -ne "s/^#% \? \?//p" Diffusion.jl > ../doc/diffusion.rst
	
../doc/schauder.rst: Schauder.jl
	sed -ne "s/^#% \? \?//p" Schauder.jl > ../doc/schauder.rst
	
../doc/nonparbayes.rst: NonparBayes.jl
	sed -ne "s/^#% \? \?//p" NonparBayes.jl > ../doc/nonparbayes.rst

../doc/randm.rst: Randm.jl
	sed -ne "s/^#% \? \?//p" Randm.jl > ../doc/randm.rst

../doc/lyap.rst: Lyap.jl
	sed -ne "s/^#% \? \?//p" Lyap.jl > ../doc/lyap.rst

../doc/linproc.rst: LinProc.jl
	sed -ne "s/^#% \? \?//p" LinProc.jl > ../doc/linproc.rst

