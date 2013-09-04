Julia package SDE.jl 
====================

This is work in progress. This package includes functionality to

* simulate diffusion processes in one or more dimension
* especially simulate vector linear processes / Ornstein-Uhlenbeck processes
* Monte Carlo sample diffusion bridges, diffusion processes conditioned to hit a point ``v`` at a prescribed time ``T``
* functions for transition density, mean and covariance of linear processes
* perform Monte Carlo estimates of transition densities of general diffusion processes
* to nonparametrically estimate the drift of a diffusion with unit diffusion coefficient

Not everything is implemented fully, the interface is crude, but most workhorse functions are there.

Layout
------

The layout/api is still the state of flux. Currently the package contains the following modules:

- **Diffusion**            Generate Ito processes and diffusions
- **SdeNonparBayes**       Nonparametrically estimate the drift of a diffusion
- **Schauder**             Provides rudimentary finite element methods and Schauder basis for SdeNonparBayes
- **Lyap**                 Computes the solution of the continuous Lyapunov equation, useful for the generation of linear processes
- **Randm**                Random symmetric, positive definite, stable matrix for testing purposes.
- **LinProc**              Homogeneous vector linear processes with additive noise


Diffusion
------------

The module `Diffusion ` contains functions to simulate Wiener processes, Wiener differentials, 
1-dimensional diffusion processes.

Worthy of mentioning are the function ``bb`` sampling Brownian bridges and the function `aug` to subsample a given Wiener process/Wiener differential 
to a higher resolution. These functions were at least tested for 1-off-errors so the Wiener processes should be *exact* 
discrete subsamples of the continuous process and not approximations valid for small time steps.


SdeNonparBayes
--------------

The producure in `SdeNonparBayes` is a Julia implementation of nonparametric Bayesian inference for
"continuously" observed one dimensional diffusion processes with unit diffusion coefficient. The drift 
is modeled as linear combination of hierarchical Faber--Schauder basis functions with a Gaussian prior 
on the coefficients. This incorporates a Brownian motion like prior on the drift function. The posterior is
then computed using Gaussian conjugacy.



Data structures
---------------

I did not introduce type definitions for stochastic processes and use vectors/arrays, so it should be easy do wrap Dataframes around everything. For the meanwhile, I like the natural notation obtained by having just vectors/arrays for dW and dt

```
N = 100
t = linspace(0., 1., N)
dt = diff(t)
X = ito(2dt + 2dW1(dt))
```

More information
----------------

[See the documentation at sdejl.readthedocs.org.](https://sdejl.readthedocs.org)

