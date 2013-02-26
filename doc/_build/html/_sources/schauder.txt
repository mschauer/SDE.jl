.. currentmodule:: Schauder
  
.. _modschauder:

Introduction to module Schauder
-------------------------------
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
---------


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
