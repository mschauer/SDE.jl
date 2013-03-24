require("SDE")
using Schauder
using Test

## procudes the same as computing mu' and then performing pick_up(mu).
function pickedup_mu(y, L)
	n = 2^(L-1) #number of even element/elements in lowest level!
	N = length(y)
	mu = zeros(2n-1)
	
	dy = [diff(y),0] #h_t = f(y_t)(x_t+h - x_t) 
	m = zeros(2n-1)
	fn = float(n)
	for t in 1:N
		for i in 1:n - 1
			fi = float(i)
			m[2i-1] = hat(y[t]*fn-fi + 1.)
			m[2i] = hat(y[t]*fn-fi +.5)
		end
		m[2n-1] = hat(y[t]*fn-fn+1)
		m = pickup_mu!(m)
		mu = mu .+ m .* dy[t]
	end
	mu
end




L = 3
N = 100

n = 2^L -1
z = rand(N)
p = finger_pm(L,0)
pm = permutationmatrix(p)

th1 = rand(n)*0.1 .+ linspace(0,1.,n)
th2 = copy(th1)
pickup!(th2)


mu1 = fe_mu(z,L, 0)
mu1b = fe_mu_c(z,L, 0)
mu2 = pickedup_mu(z, L)
mu2b = pickup_mu!(copy(mu1))

tic()
Sigma1 = fe_Sigma_at(z,0.1, L)
toc()
tic()
Sigma1b = fe_Sigma_dot(z,0.1, L)
toc()

Sigma2 = copy(Sigma1)
Sigma2 = pickup_Sigma!(Sigma2)

println([th1 th2 mu1 mu1b mu2 mu2b])

println("Should be zero:")
println(norm(mu2-mu2b))
println(norm(mu1-mu1b))
println(norm(Sigma1-Sigma1b))
println("Should be equal:")

println(th1'*mu1, th2'*mu2, th2'*mu2b)

println(th1'*Sigma1*th1,th1'*Sigma1b*th1, th2'*Sigma2*th2)

