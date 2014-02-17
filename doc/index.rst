.. SDE documentation master file, created by
   sphinx-quickstart on Sun Feb 24 21:09:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Julia package SDE.jl 
====================

Layout
------

The main module 

*SDE*
    Simulate diffusion processes in one or more dimension.
    Especially simulate vector linear processes / Ornstein-Uhlenbeck processes
    Monte Carlo sample diffusion bridges, diffusion processes conditioned to hit a point v at a prescribed time T
    Functions for transition density, mean and covariance of linear processes
    Perform Monte Carlo estimates of transition densities of general diffusion processes

with a Submodule

*SDE.Schauder*
    To nonparametrically estimate the drift of a diffusion with unit diffusion coefficient
    using a Schauder wavelet "basis"

The package contains the additional modules:

*Diffusion*
   Alternative API to generate Ito processes and diffusions

*Randm*
   Random symmetric, positive definite, stable matrix for testing purposes.


Method
------

The producure in ``SDE.Schauder`` is a Julia implementation of nonparametric Bayesian inference for
"continuously" observed one dimensional diffusion processes with unit diffusion coefficient. The drift 
is modeled as linear combination of hierarchical Faber--Schauder basis functions with a Gaussian prior 
on the coefficients. This incorporates a Brownian motion like prior on the drift function. The posterior is
then computed using Gaussian conjugacy.

This is work in progress.


Location of the documentation
-----------------------------

https://sdejl.readthedocs.org


Contents:

.. toctree::
   :maxdepth: 2

   SDE   
   schauder
   diffusion
   randm


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

