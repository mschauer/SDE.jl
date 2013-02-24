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
