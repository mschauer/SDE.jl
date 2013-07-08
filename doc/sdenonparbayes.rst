.. currentmodule:: SdeNonparBayes
  
Introduction to module SdeNonparBayes
-------------------------------------

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


Optional additial basis functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
---------

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
.. function:: visualize_posterior(post[, truedrift])
	
	Plot 2r*se wide marginal credibility bands, where ``post`` is the result of 
	bayes_drift and truedrift the true drift (if known :-) ).

