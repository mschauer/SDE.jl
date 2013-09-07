using SDE
using Lyap
using LinProc
using Winston
using Randm
using Optim

srand(7)
d = 2 # dimension
M = 500 #number of observations
TT = .1*M #total time span
function norma!(x)
 minx = min(x)
 x = (x - minx) / (max(x)-minx)
 x
end
dt = diff(TT*norma!(sort!(rand(M))))


# diffusion coefficient
sigma = [ 0.9   .2; -.2  .7] 
A = sigma*sigma'

# starting point
u = [0.,0.] 

# true drift b(x) = B0*x + beta0
B0 = [  -0.2  -1.0;  
	 0.5  -0.4]
beta0 = [0.,0.] 

# many functions in LinProc expect the solution to the lyapunov equation given B0 and A as argument instead of A
lambda0 = lyap(B0', -A)

# simulate process
function linexact(u, B, beta, lambda, dt)
	X = zeros(d, M)
	X[:,1] = u
	for i in 1 : M-1
		# sample from the transition probability
		X[:,i+1] = LinProc.sample_p(dt[i], X[:,i], B, beta, lambda) 
	end
	X
end

# compute log likelihood
function llikeli(X, B, beta, lambda, dt)
	ll = 0.0
	for i in 1 : M-1
		ll += LinProc.lp(dt[i], X[:,i], X[:,i+1], B, beta, lambda) 
	end
	ll
end



# simulate exact Ornstein--Uhlenbeck process with parameter B0, beta0, A (lambda0)
X =  linexact(u, B0, beta0, lambda0, dt)

# likelihood of true B
llB0 =  (llikeli(X, B0, beta0, lambda0, dt))
println("true B0")
display(B0)
println("\nlikelihood of true B $llB0")

# looks like this
plot(X[1,:], X[2,:])

# stable d-dim matrix parametrized with 2d^2 numbers, real eigenvalues < 0.01
function stable(Y, d=2)
	# positive definite matrix
	x = reshape(Y[1:d*d], 2,2)
	a = x'*x 
	
	# anti symmetric matrix
	above_diag =	[ i<j ? 1 : 0 for i in 1:d, j in 1:d] # nonzero pattern for above diagonal
	y = reshape(Y[d*d+1 : 2*d*d], 2,2)
	b = above_diag.*y - above_diag'.*y' 
	
	# return stable matrix
	b - a -0.01eye(2)
end

function objective(Y)
	assert(size(Y) == (8,))
	B = stable(Y) # obtain stable matrix corresponding to numbers Y
	lambda = lyap(B', -A)	# lambda depends on B
	ll = try
		-llikeli(X, B, beta0, lambda, dt)
	catch y
		println(y)
		Inf
	end
	ll
end


#B=round(ml(500),3) 
O = optimize(objective, ones(8))
print(O)
B = stable(O.minimum)
llmax = -O.f_minimum
println(llmax, " ", llB0)
println([B B0])
