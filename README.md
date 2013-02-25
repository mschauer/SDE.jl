Itostat package
===============

Layout
------
The package contains the following modules:

- **Diffusion**            Generate Ito processes and diffusions
- **SdeNonparBayes**       Nonparametrically estimate the drift of a diffusion
- **Schauder**             Provides rudimentary finite element methods and Schauder basis for SdeNonparBayes

Method
------

The producure in `SdeNonparBayes` is a Julia implementation of nonparametric Bayesian inference for
"continuously" observed one dimensional diffusion processes with unit diffusion coefficient. The drift 
is modeled as linear combination of hierarchical Faber--Schauder basis functions with a Gaussian prior 
on the coefficients. This incorporates a Brownian motion like prior on the drift function. The posterior is
then computed using Gaussian conjugacy.

This is work in progress.



Packaging steps remaining
-------------------------
- From within Julia, execute the command `Pkg.version("Packagename", v"0.0.0")`. This will  create a directory `.julia/METADATA/Packagename`, which is what you need to have your package registered.
- In `.julia/METADATA/Packagename`, create a file `url` that contains the *read-only* URL for your public repository. For example, git://github.com/mygithubaccount/Packagename.jl.git.
- Commit your changes to `.julia/METADATA`, push to your public `METADATA` repository, and submit a pull request.




