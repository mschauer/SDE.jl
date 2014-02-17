all: ../doc/diffusion.rst  ../doc/randm.rst ../doc/SDE.rst ../doc/schauder.rst

../doc/diffusion.rst: Diffusion.jl
	sed -ne "s/^#% \? \?//p" Diffusion.jl > ../doc/diffusion.rst
	
../doc/randm.rst: Randm.jl
	sed -ne "s/^#% \? \?//p" Randm.jl > ../doc/randm.rst

../doc/SDE.rst: SDE.jl
	sed -ne "s/^#% \? \?//p" SDE.jl > ../doc/SDE.rst

../doc/schauder.rst: Schauder.jl npbayes.jl
	sed -ne "s/^#% \? \?//p" Schauder.jl > ../doc/schauder.rst
	sed -ne "s/^#% \? \?//p" npbayes.jl >> ../doc/schauder.rst


