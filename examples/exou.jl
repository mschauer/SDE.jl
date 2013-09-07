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
pritnln(B0)
println("likelihood of true B ", round(llB0,3))

# looks like this
plot(X[1,:], X[2,:])

# we need to parametrize stable matrices
# matrices with eigenvalues with strictly negative real parts

# stable d-dim matrix parametrized with 2d^2 numbers, real eigenvalues < 0.01
function stable(Y, d=2, ep=0.01)

	# convert first d*(d+1)/2 values of Y into upper triangular matrix
	# positive definite matrix
	x = zeros(d,d)
	k = 1
	for i in 1:d
		for j in i:d
		x[i,j] = Y[k]
		k = k + 1
		end
	end
	# convert next d*(d+1)/2 -d values of Y into anti symmetric matrix
	y = zeros(d,d)
	for i in 1:d
		for j  in i+1:d
		y[i,j] = Y[k]
		y[j,i] = -y[i, j]
		k = k + 1
		end
	end
	assert(k -1 == d*d == length(Y))
	
	# return stable matrix as a sum of a antisymmetric and a positive definite matrix
	y - x'*x - ep*eye(2) 
end

function objective(Y)
	assert(length(Y) == d*d)
	B = stable(Y) # obtain stable matrix corresponding to numbers Y
	lambda = lyap(B', -A)	# lambda depends on B
	ll = try
		-llikeli(X, B, beta0, lambda, dt)
	catch y
		if isa(y, Base.PosDefException)
			println("Skip numerical indefinite matrix .")
			1E16			
		elseif isa(y, Base.Singular)
			println("Skip numerical singular matrix.")
			1E16			
		else 
		 throw(y)
		end

	end
	ll
end


#B=round(ml(500),3) 
O = optimize(objective, ones(4))
print(O)
B = round(stable(O.minimum),3)
llmax = -O.f_minimum
println("Estimated B (llikelihood ",round(llmax,3),")\n",round(B,3), "\nTrue B0 (llikelihood ", round(llB0,3),")\n",B0 )
