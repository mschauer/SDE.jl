.. currentmodule:: Randm
  
Introduction to module Randm
----------------------------

Random matrices for testing purposes. I did not figure out the actual distributions
the matrices are drawn from.

Reference 
---------

.. function:: randsym(d)
             
	Random symmetric matrix.
	
.. function:: randposdef(d)
             
	Random positive definite matrix of dimension ``d``.
	
.. function:: randstable(d)
             
	Random stable matrix (matrix with eigenvalues with negative real part) with
	dimension ``d``.
