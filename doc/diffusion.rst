.. currentmodule:: Diffusion
  
.. |Ito| unicode:: Ito

Module Diffusion
----------------

Introduction
~~~~~~~~~~~~

The functions in this module operate on three conceptual different objects, (although they 
are currently just represented as vectors and arrays.)

Stochastic processes, denoted ``x``, ``y``, ``w`` are arrays of values which 
are sampled at distance ``dt``, ``ds``, where ``dt``, ``ds`` are either scalar or 
Vectors ``length(dt)=size(W)[end]``.
Stochastic differentials are denoted ``dx``, ``dw`` etc., and are first differences of 
stochastic processes. Finally, t can denote the total time or correspond to a vector of
``size(W)[end]`` sampling time poins.

Note the following convention: In analogy with the definition of the |Ito| integral,

	intxdw[i] = x[i]](w[i+1]-w[i]) (== x[i]dw[i])

and

	length(w) = length(dw) + 1

Reference 
~~~~~~~~~

.. function:: brown1(u, t, n::Integer) 

	Compute ``n`` equally spaced samples of 1d Brownian motion in
	the interval ``[0,t]``, starting from point ``u``

.. function:: brown(u, t, d::Integer, n::Integer) 

	Simulate ``n`` equally spaced samples of ``d``-dimensional Brownian motion in
	the interval ``[0,t]``, starting from point ``u``

.. function:: dW1(t, n::Integer)
              dW(t, d::Integer, n::Integer)

	Simulate a ``1``-dimensional (``d``-dimensional)
	Wiener differential with ``n`` values in the 
	the interval ``[0,t]``, starting from point ``u``
	
.. function:: dW(dt::Vector, d::Integer) 

	Simulate a ``d``-dimensional Wiener differential sampled at
	time points with distances given by the vector ``dt``		
	
.. function:: ito(y, dx)
              ito(dx)
              cumsum0(dx)

	Integrate a valued stochastic process with respect to a stochastic differential.
	R, R^2 (d rows, n columns), R^3.
	
	``ito(dx)`` is a shortcut for ``ito(ones(size(dx)[end], dx)``.
	So ``ito(dx)`` is just a ``cumsum0`` function which is a inverse to ``dx = diff([0, x1, x2, x3,...])``.

.. function:: ..(y, dx)
              ydx(y, dx)

	``y .. dx`` returns the stochastic differential ``ydx`` defined by the property

		ito(ydx) == ito(y, dx)

.. function:: bb(u, v, t, n) 

	Simulates ``n`` equidistant samples of a Brownian bridge from point ``u`` to ``v`` in time ``t``

.. function:: dWcond1(v,t,n)

  	Simulates ``n`` equidistant samples of a "bridge noise": that is a Wiener differential ``dW``
  	conditioned on ``W(t) = v``


.. function:: aug(dw, dt, n)
              aug(dt, n)

	Take Wiener differential sampled at ``dt`` and return Wiener differential subsampled ``n`` times
	between each observation with new length ``length(dw)*n``.
	``aug(dt,n)`` computes the corresponding subsample of times.

.. function:: quvar(x)
             
	Computes quadratic variation of ``x``.
	
.. function:: bracket(x)
              bracket(x,y)

	Computes quadratic variation process of ``x`` (of ``x`` and ``y``).
	
.. function:: euler(t0, u, b, sigma, dt, dw)
              euler(t0, u, b, sigma, dt)

	Simulates a 1-dimensional diffusion process using the Euler-Maruyama approximation
	with drift ``b(t,x)`` and diffusion coefficient ``sigma(t,x)``
	starting in ``(t0, u)`` using ``dt`` and given Wiener differential ``dw``.

