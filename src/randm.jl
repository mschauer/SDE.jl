#some random matrices for testing
module Randm
export randsym, randposdef, randstable, randnormal, randorth, randunitary

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

##%  .. function:: randsym(d)
##%               
##%  	Random symmetric matrix. Note: randsym
##%  	

#function randsym(d)
#	x = 0.5*randn(d, d)
#	x' + x
#end

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

#%  .. function:: randunitary(d)
#%               
#%  	Random unitary matrix of dimension ``d``.
#%  	

randunitary(d) = expm(im* randposdef(d))

#%  .. function:: randorth(d)
#%               
#%  	Orthogonal matrix drawn according to the Haar measure on the group of orthogonal matrices.
#%  	

randorth(d) = qr(randn(d,d))[1]

#%  .. function:: randnormal(d) 
#%               
#%  	Random normal matrix.
#%  	


function tridiag(di, u, l)
	assert(length(di) -1 == length(u) == length(l))
	M = diagm(di)
	for i in 1:length(u)
	  M[i,i+1] = l[i]
  	  M[i+1,i] = u[i]
	end
	M
	
end
tridiag(alpha, beta) = tridiag(alpha, beta, -beta)


function randnormal(d) 
	Q = randorth(d)
	B = schur(randn(d,d))[1] #same eigenvalues as a randn(d,d)
	m = d
	alpha = zeros(d)
	beta = zeros(d-1)
	while m > 1
	   s = abs(B[m-1,m-1]) + abs(B[m,m])
	   if s + abs(B[m,m-1]) > s # if significant offdiagonal value: evaluate submatrix m-1:m;
		spur = B[m,m] + B[m-1,m-1]
		dis2 = square(B[m,m]- B[m-1,m-1]) + 4.0*B[m,m-1]*B[m-1,m] #= spur^2 - 4det
		alpha[m] = alpha[m-1] = 0.5*spur
		
		beta[m-1] = 0.5*sqrt(-dis2)
		m -= 2       
	   else
	        alpha[m] = B[m,m]
	        beta[m-1] = 0
	        m -= 1        
	   end
	   
	end
	if (m == 1) 
		alpha[1] = B[1,1]
	end
	
	Q' * tridiag(alpha, beta) * Q
	#tridiag(alpha, beta)
	
end

end #module
