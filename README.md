Julia package SDE.jl 
====================

This is work in progress. This package includes functionality to

* simulate diffusion processes in one or more dimension
* especially simulate vector linear processes / Ornstein-Uhlenbeck processes
* Monte Carlo sample diffusion bridges, diffusion processes conditioned to hit a point *v* at a prescribed time *T*
* functions for transition density, mean and covariance of linear processes
* perform Monte Carlo estimates of transition densities of general diffusion processes
* to nonparametrically estimate the drift of a diffusion with unit diffusion coefficient

The layout/api is still the state of flux and the package is currently undergoing a major revision,
I hope to have some nice version 1 in summer.

[![Build Status](https://api.travis-ci.org/mschauer/SDE.jl.png)](https://travis-ci.org/mschauer/SDE.jl)