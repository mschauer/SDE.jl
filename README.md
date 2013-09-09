Julia package SDE.jl 
====================

This is work in progress. This package includes functionality to

* simulate diffusion processes in one or more dimension
* especially simulate vector linear processes / Ornstein-Uhlenbeck processes
* Monte Carlo sample diffusion bridges, diffusion processes conditioned to hit a point *v* at a prescribed time *T*
* functions for transition density, mean and covariance of linear processes
* perform Monte Carlo estimates of transition densities of general diffusion processes
* to nonparametrically estimate the drift of a diffusion with unit diffusion coefficient

Not everything is implemented fully, the interface is crude, but most workhorse functions are there.

Overview
--------

The layout/api is still the state of flux. Currently the package contains the following modules:

- **Diffusion**            Generate Ito processes and diffusions
- **NonparBayes**          Nonparametrically estimate the drift of a diffusion
- **Schauder**             Provides rudimentary finite element methods and Schauder basis for NonparBayes
- **Lyap**                 Computes the solution of the continuous Lyapunov equation, useful for the generation of linear processes
- **Randm**                Random symmetric, positive definite, stable matrix for testing purposes.
- **LinProc**              Homogeneous vector linear processes with additive noise

[See the documentation at sdejl.readthedocs.org.](https://sdejl.readthedocs.org)



Module Diffusion
----------------

The module `Diffusion ` contains functions to simulate Wiener processes, Wiener differentials, 
1-dimensional diffusion processes.

Worthy of mentioning are the function ``bb`` sampling Brownian bridges and the function `aug` to subsample a given Wiener process/Wiener differential 
to a higher resolution. These functions were at least tested for 1-off-errors so the Wiener processes should be *exact* 
discrete subsamples of the continuous process and not approximations valid for small time steps.


Module NonparBayes
------------------

The module  `NonparBayes` is a Julia implementation of nonparametric Bayesian inference for
"continuously" observed one dimensional diffusion processes with unit diffusion coefficient. The drift 
is modeled as linear combination of hierarchical Faber--Schauder basis functions with a Gaussian prior 
on the coefficients. This incorporates a Brownian motion like prior on the drift function. The posterior is
then computed using Gaussian conjugacy.



Example for LinProc
-------------------

The following example gives a look and feel shows how the Module `LinProc` could be used. This is taken from `examples/exou.jl`.

We first simulate a two dimensional Ornstein--Uhlenbeck process.

```julia
using SDE
using Lyap
using LinProc

using Winston

srand(7)
d = 2 # dimension
M = 500 #number of observations
TT = .1*M #total time span

# normalize values, returns values in interval [0,1]
function norma!(x)
 minx = min(x)
 x = (x - minx) / (max(x)-minx)
 x
end

grid = TT*norma!(sort!(rand(M))) # a grid of random design points in the interval [0, TT]

# diffusion coefficient (assumed to be known in this example)
sigma = [ 0.9   .2; -.2  .7] 
A = sigma*sigma'

# starting point
u = [0.,0.] 

# true drift function, b(x) = B0*x + beta0 (unknown, used to obtain observations)

B0 = [  -0.2  -1.0;  # true mean reversion matrix
	 0.5  -0.4]
beta0 = [0.,0.] 

# many functions in LinProc expect the solution to the 
# Lyapunov equation given B and A as argument instead of A 

lambda0 = lyap(B0', -A) # solves B*lambda + lambda*B' = -A

# and many functions need the distances between the design points

dt = diff(grid) 


# simulate exact Ornstein--Uhlenbeck process with parameter B0, beta0, A (lambda0)
X =  linexact(u, B0, beta0, lambda0, dt)

# looks like this
println("Plotting observations")
plot(X[1,:], X[2,:])

```
![Julia plot of X](https://raw.github.com/mschauer/SDE.jl/master/doc/exou.jl.png)

We now proceed to estimate the mean reversion matrix `B0` from the observations:

```julia
using Optim

# provide objective function for maximum likelihood

function objective(Y)

	assert(length(Y) == d*d) # Y is a vector of length d*d parametrizing all stable matrices B	
	B = LinProc.stable(Y, d, 0.02) # obtain stable matrix with eigenvalues with real part < 0.02 corresponding to numbers in Y
	lambda = lyap(B', -A)	# lambda depends on B

	ob =  #minimize likelihood
	try
		-linll(X, B, beta0, lambda, dt) # negative discrete observations loglikelihood
	catch y #catch numerical singularies
		if isa(y, Base.PosDefException)
			println("Skip numerical indefinite matrix")
			1.0E16 # move away! 			
		elseif isa(y, Base.Singular)
			println("Skip numerical singular matrix")
			1.0E16	# move away!		
		else 
		 throw(y)
		end

	end
	ob
end

# find maximum likelihood estimate for B
println("Maximize likelihood...") 
O = optimize(objective, ones(4))
println(O)
B = round(LinProc.stable(O.minimum,2,0.02),3)
llmax = -O.f_minimum
# as comparison: loglikelihood of true B
llB0 =  linll(X, B0, beta0, lambda0, dt)

print("\nEstimated mean reversion matrix B (log-likelihood ",round(llmax,3),")\n",round(B,3), "\nTrue B0 (log-likelihood ", round(llB0,3),")\n",B0 )
```

Running the program gives
```
Estimated mean reversion matrix B (log-likelihood 236.051)
-.216	-.973
.431	-.394

True B0 (log-likelihood 235.486)
-.2	-1
.5	-.4
```

More information
----------------

[See the documentation at sdejl.readthedocs.org.](https://sdejl.readthedocs.org)

[![Build Status](https://api.travis-ci.org/mschauer/SDE.jl.png)](https://travis-ci.org/mschauer/SDE.jl)
