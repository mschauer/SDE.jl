# Example: Maximum-likelihood estimate of the 

using SDE
using Lyap
using LinProc

using Winston
using Optim

srand(7)
d = 2 # dimension
M = 500 #number of observations
TT = .1*M #total time span

# normalize values, returns values in interval [0,1]
function norma!(x)
 minx = min(x)
 x = (x - minx) / (max(x)-minx)
 x
end

grid = TT*norma!(sort!(rand(M))) # a grid of random design points in the interval [0, TT]

# diffusion coefficient (assumed to be known in this example)
sigma = [ 0.9   .2; -.2  .7] 
A = sigma*sigma'

# starting point
u = [0.,0.] 

# true drift function, b(x) = B0*x + beta0 (unknown, used to obtain observations)

B0 = [  -0.2  -1.0;  # true mean reversion matrix
	 0.5  -0.4]
beta0 = [0.,0.] 

# many functions in LinProc expect the solution to the 
# Lyapunov equation given B and A as argument instead of A 

lambda0 = lyap(B0', -A) # solves B*lambda + lambda*B' = -A

# and many functions need the distances between the design points

dt = diff(grid) 


# simulate exact Ornstein--Uhlenbeck process with parameter B0, beta0, A (lambda0)
X =  linexact(u, B0, beta0, lambda0, dt)

# as comparison: loglikelihood of true B
llB0 =  linll(X, B0, beta0, lambda0, dt)
println("True mean reversion matrix B0")
print(B0)
println("Likelihood of B0 ", round(llB0,3))

# looks like this
println("Plotting observations")
plot(X[1,:], X[2,:])


# objective function for maximum likelihood

function objective(Y)

	assert(length(Y) == d*d) # Y is a vector of length d*d parametrizing all stable matrices B	
	B = LinProc.stable(Y, d, 0.02) # obtain stable matrix with eigenvalues with real part < 0.02 corresponding to numbers in Y
	lambda = lyap(B', -A)	# lambda depends on B

	ob =  #minimize likelihood
	try
		-linll(X, B, beta0, lambda, dt) # negative discrete observations loglikelihood
	catch y #catch numerical singularies
		if isa(y, Base.PosDefException)
			println("Skip numerical indefinite matrix")
			1.0E16 # move away! 			
		elseif isa(y, Base.Singular)
			println("Skip numerical singular matrix")
			1.0E16	# move away!		
		else 
		 throw(y)
		end

	end
	ob
end

# find maximum likelihood estimate for B
println("Maximize likelihood...") 
O = optimize(objective, ones(4))
println(O)
B = round(LinProc.stable(O.minimum,2,0.02),3)
llmax = -O.f_minimum
print("\nEstimated mean reversion matrix B (log-likelihood ",round(llmax,3),")\n",round(B,3), "\nTrue B0 (log-likelihood ", round(llB0,3),")\n",B0 )
