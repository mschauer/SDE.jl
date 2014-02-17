.. currentmodule:: Schauder
  
.. _modschauder:

Module Schauder
---------------

Introduction
~~~~~~~~~~~~

In the following ``hat(x)`` is the piecewise linear function taking values
values ``(0,0), (0.5,1), (1,0)`` on the interval ``[0,1]`` and ``0`` elsewhere.

The Schauder basis of level ``L > 0`` in the interval ``[a,b]`` can be defined recursively
from ``n = 2^L-1`` classical finite elements :math:`\psi_i(x)` on the grid

	``a + (1:n)/(n+1)*(b-a)``.


Assume that f is expressed as linear combination 

	:math:`f(x) = \sum_{i =1}^n c_i \psi_i(x)`

with 

	:math:`\psi_{2j-1}(x) = hat(nx-j + 1)` 	for :math:`j = 1 \dots 2^{L-1}` 

and 

	:math:`\psi_{2j}(x) = hat(nx-j + 1/2)` 	for :math:`j = 1 \dots 2^{L-1}-1`

Note that these coefficients are easy to find for the finite element basis, just

.. code-block:: matlab

	function fe_transf(f, a,b, L)
	    n = 2^L-1 	
	    return map(f, a + (1:n)/(n+1)*(b-a))
	end

  
Then the coefficients of the same function with respect to the Schauder basis 

	:math:`f(x) = \sum_{i = 1}^n c_i \phi_i(x)`

where for ``L = 2``

	:math:`\phi_2(x) = 2 hat(x)`

	:math:`\phi_1(x) = \psi_1(x) = hat(2x)`

	:math:`\phi_3(x) = \psi_3(x) = hat(2x-1)`

can be computed directly, but also using the recursion

	:math:`\phi_2(x) = \psi_1(x) + 2\psi_2(x) + \psi_2(x)`

This can be implemented inplace (see ``pickup()``), and is used throughout. 


.. code-block:: matlab

	for l in 1:L-1
		b = sub(c, 2^l:2^l:n)
		b[:] *= 0.5
		a = sub(c, 2^(l-1):2^l:n-2^(l-1))
		a[:] -= b[:]
		a = sub(c, 2^l+2^(l-1):2^l:n)
		a[:] -= b[1:length(a)]
	end

Reference 
~~~~~~~~~


.. function:: pickup!(x)

	Inplace computation of 2^L-1 Schauder-Faber coefficients from 
	``2^L-1`` overlapping finite-element coefficients ``x``.

	-- inverse of ``Schauder.drop``

	-- L = level(xj)

.. function:: drop!(x)

	Inplace computation of 2^L-1 finite element coefficients from 
	2^L-1 Faber schauder coefficients ``x``.

	-- inverse of ``Schauder.pickup``

.. function:: finger_permute(x)

	Reorders vector ``x`` or matrix ``A`` according to the reordering
	of the elements of a Faber-Schauder-basis from
	left to right, from bottom to top.

.. function:: finger_pm(L, K)

	Returns the permuation used in ``finger_permute``.
	Memoized reordering of faber schauder elements from low level to high level. The last K elements/rows are left untouched.

.. function:: level(x)

	Gives the no. of levels of the biggest Schauder basis with less then length(x) elements.
		level(x) = ilogb(size(x,1)+1)

.. function:: level(x, K)

	Gives the no. of levels ``l`` of the biggest Schauder basis with less then length(x) elements
	and the number of additional elements ``n-2^l+1``.

.. function:: vectoroflevels(L, K)
             
	Gives a vector with the level of the hierarchical elements.
	
.. function:: hat(x) 
             
	Hat function. Piecewise linear functions with values (-inf,0), (0,0),(0.5,1), (1,0), (inf,0).
	-- ``x`` vector or number
Introduction
~~~~~~~~~~~~


The procedure is as follows.
Consider the diffusion process :math:`(x_t\colon 0 \le t \le T)` given by

	:math:`dx_t = b(x_t) dt + dw_t`


where the drift ``b`` is expressed as linear combination 

	:math:`f(x) = \sum_{i =1}^n c_i \phi_i(x)`

(see :ref:`modschauder`) and 
prior distribution on the coefficients 

	:math:`c_i \sim N(0,\xi_i)`

Then the posterior distribution of :math:`b` given observations :math:`x_t` is given by

	:math:`c_i | x_s \sim N(W^{-1}\mu, W^{-1})`
	:math:`W = \Sigma + (\operatorname{diag}(\xi))^{-1}`, 

with the nxn-matrix 

	:math:`\Sigma_{ij} = \int_0^T \phi_i(x_t)\phi_j(x_t) dt`

and the n-vector

	:math:`\mu_i = \int_0^T \phi_i(x_t) d x_t`.


Using the recursion detailed in :ref:`modschauder`, one rather computes

	:math:`\Sigma^\prime_{ij} = \int_0^T \psi_i(x_t)\psi_j(x_t) dt`

and the n-vector

	:math:`\mu^\prime_i = \int_0^T \psi_i(x_t) d x_t`

and uses ``pickup_mu!(mu)`` and ``pickup_Sigma!(Sigma)`` to obtain :math:`\mu` and :math:`\Sigma`. 


Optional additional basis functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can extend the basis by additional functions, implemented are variants. ``B1`` includes a constant, ``B2`` two linear functions

``B1``
	:math:`\phi_1 \dots \phi_n, c`

``B2``
	:math:`\phi_1 \dots \phi_n, \max(1-x, 0), \max(x, 0)`

To compute ``mu``, use

.. code-block:: matlab

	mu = pickup_mu!(fe_mu(y,L, 0))
	mu = fe_muB1(mu, y);

or

.. code-block:: matlab

	mu = pickup_mu!(fe_mu(y,L, 0))
	mu = fe_muB2(mu, y);


Reference 
~~~~~~~~~

Functions taking ``y` without parameter [a,b] expect ``y`` to be shifted into the intervall ``[0,1]``.

.. function:: pickup_mu!(mu)
             
	computes mu from mu'
	
.. function:: drop_mu!(mu)
             
	Computes mu' from mu.
	
.. function:: pickup_Sigma!(Sigma)
             
	Transforms Sigma' into Sigma.
	
.. function:: drop_Sigma!(Sigma)
             
	Transforms Sigma into Sigma'.
	
.. function:: fe_mu(y, L, K)
             
	Computes mu' from the observations `y` using ``2^L-1`` basis elements
	and returns a vector with ``K`` trailing zeros (in case one ones to customize
	the basis.
.. function:: fe_muB1(mu, y)
             
	Append :math:`\mu_{n+1} = \int_0^T \phi_{n+1} d x_t` with :math:`\phi_{n+1} = 1`.
	
.. function:: fe_muB2(mu, y)
             
	Append :math:`\mu_{n+1} = \int_0^T \phi_{n+1} d x_t` with :math:`\phi_{n+1} =  \max(1-x, 0)`
	and :math:`\mu_{n+2} = \int_0^T \phi_{n+2} d x_t` with :math:`\phi_{n+2} =  \max(x, 0)`

.. function:: fe_Sigma(y, dt, L)
             
	Computes the matrix Sigma' from the observations `y` uniformly spaced at distance ``dt``
	using ``2^L-1`` basis elements.

.. function:: bayes_drift(x, dt, a, b, L, xirem, beta, B)
             
	Performs estimation of drift on observations ``x`` in [a,b] spaced at distance ``dt``
	using the Schauder basis of level ``L`` and level wise coefficients decaying at rate ``beta``.
	A Brownian motion like prior is obtained for beta= 0.5. The ``K`` remaining optional 
	basiselements have variance ``xirem``. 

	The result is returned as ``[designp coeff se]`` where ``coeff`` are coefficients of finite elements with maximum at the designpoints ``designp`` and standard error ``se``.


	
	Observations outside [a,b] may influence the result through ``phi_{n+1}, ..., phi_{n+K}``
