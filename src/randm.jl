#some random matrices for testing
module Randm
export randsym, randposdef, randstable

#%  .. currentmodule:: Randm
#%    

#%  Introduction to module Randm
#%  ----------------------------
#%  
#%  Random matrices for testing purposes. I did not figure out the actual distributions
#%  the matrices are drawn from.
#%  


#%  Reference 
#%  ---------
#%  

#%  .. function:: randsym(d)
#%               
#%  	Random symmetric matrix.
#%  	

function randsym(d)
	x = randn(d, d)
	x' + x
end

#%  .. function:: randposdef(d)
#%               
#%  	Random positive definite matrix of dimension ``d``.
#%  	

function randposdef(d)
	x = randn(d, d) 
	x'* x
end

#%  .. function:: randstable(d)
#%               
#%  	Random stable matrix (matrix with eigenvalues with negative real part) with
#%  	dimension ``d``.


function randstable(d)
	# positive definite matrix
	x = randn(d, d) 
	a = x'*x 
	
	# anti symmetric matrix
	above_diag =	[ i<j ? 1 : 0 for i in 1:d, j in 1:d] # nonzero pattern for above diagonal
	y = randn(d, d) 
	b = above_diag.*y - above_diag'.*y' 
	
	# return stable matrix
	b - a
end



end
