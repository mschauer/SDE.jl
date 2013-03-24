Julia package SDE.jl 
====================

Layout
------
The package contains the following modules:

- **Diffusion**            Generate Ito processes and diffusions
- **SdeNonparBayes**       Nonparametrically estimate the drift of a diffusion
- **Schauder**             Provides rudimentary finite element methods and Schauder basis for SdeNonparBayes
- **Lyap**                 Computes the solution of the continuous Lyapunov equation, useful for the generation of linear processes
- **Randm**                Random symmetric, positive definite, stable matrix for testing purposes.
- **LinProc**              Homogeneous vector linear processes with additive noise


Method
------

The producure in `SdeNonparBayes` is a Julia implementation of nonparametric Bayesian inference for
"continuously" observed one dimensional diffusion processes with unit diffusion coefficient. The drift 
is modeled as linear combination of hierarchical Faber--Schauder basis functions with a Gaussian prior 
on the coefficients. This incorporates a Brownian motion like prior on the drift function. The posterior is
then computed using Gaussian conjugacy.

This is work in progress.


Data structures
---------------

I did not introduce type definitions for stochastic processes and use vectors/arrays, so it should be easy do wrap Dataframes around everything. For the meanwhile, I like the natural notation obtained by having just vectors/arrays for dW and dt

	N = 100
	t = linspace(0., 1., N)
	dt = diff(t)
	X = ito(2dt + 2dW1(dt))

More information
----------------

[See the documentation.](doc/_build/html/)

