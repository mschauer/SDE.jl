.. currentmodule:: Randm
  
Introduction to module Randm
----------------------------

Random matrices for testing purposes. I did not figure out the actual distributions
the matrices are drawn from.

Reference 
---------

.. function:: randposdef(d)
             
	Random positive definite matrix of dimension ``d``.
	
.. function:: randstable(d)
             
	Random stable matrix (matrix with eigenvalues with negative real part) with
	dimension ``d``.
.. function:: randunitary(d)
             
	Random unitary matrix of dimension ``d``.
	
.. function:: randorth(d)
             
	Orthogonal matrix drawn according to the Haar measure on the group of orthogonal matrices.
	
.. function:: randnormal(d) 
             
	Random normal matrix.
	
