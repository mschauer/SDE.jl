# test linproc.jl
#require("SDE")

using Cubature
#using LinProc
using Randm
#using Lyap
using Base.Test


d= 4
# B = [-10.3268 .67701   -1.85162;
#   2.60125  -5.35212    0.4388; 
#  -1.42184   0.904901  -0.5423]
#
srand(5)
B = randstable(d)
A = randposdef(d)
function aa(s)
 A
end
lambda = Lyap.lyap(B', -A)

T = 1
t = 0.5
v = zeros(d*d)
function Hquad(t, T, B, a::Function)
	function f(s, v)
		v[:] = vec(expm(-(s-t)*B)*a(s)*expm(-(s-t)*B)')
	end
	Q = reshape(hquadrature(d*d, f, t, T, 1E-15, 1E-15, 5000)[1], d,d)
	 -inv(Q)
end	

println("Test LinProc.H, LinProc.K")

H1 = LinProc.H(T-t, B, lambda)
H2 = Hquad(t, T, B, aa)
K1 = LinProc.K(T-t, B, lambda)
H3 = -expm((T-t)*B')*inv(K1)*expm((T-t)*B)


h = T-t
phi = expm(h*B)
phim = expm(-h*B)
K2 = lambda - phi*lambda*phi'
H4 = -phi'*inv(lambda - phi*lambda*phi')*phi
H5 = LinProc.H(h, B, lambda)

@test norm(H1 - H2) < 1E-10
@test norm(H1 - H3) < 1E-10
@test norm(H1 - H4) < 1E-10
@test norm(H1 - H5) < 1E-10

@test norm(K1-K2) < 16*eps()

#obtaining r via quadrature
function varmu(h, x, B, beta)
	function f(s, y)
		y[:] = expm(-s*B)*beta
	end
	integral = Cubature.hquadrature(d, f, 0, h, 1E-15, 1E-15, 800)[1]
	expm(h*B)*(x + integral)
end	


function varr(h, x, v, B, beta, lambda)
	mu = varmu(h, x, B, beta)
	expm(h*B')*inv(LinProc.K(h, B, lambda))*(v - mu) 
end

function varr3(h, x, v, B, beta, lambda)
	binv = inv(B)
	phi = expm(h*B)
	phim = expm(-h*B)
	vt =  phim*( v) + phim * binv * beta -binv*beta 
	LinProc.H(h, B, lambda) * (x - vt )
	
end

#H3 = -expm((T-t)*B')*inv(K1)*expm((T-t)*B)

println("Test LinProc.r")

srand(5)
x = randn(d)
v = randn(d)
beta = randn(d)
r1 = LinProc.r(T-t, x, v, B, beta, lambda)
r2 = varr(T-t, x, v, B, beta, lambda)
r3 = varr3(T-t, x, v, B, beta, lambda)

@test norm(r1 - r2) < 1E-10
@test norm(r1 - r3) < 1E-10

println("Test LinProc.mu")

mu1 = LinProc.mu(h, x, B, beta)
mu2 = varmu(h, x, B, beta)

@test norm(mu1 - mu2) < 1E-14


#println(mu1)
#println(mu2)

#checks if transition density p integrates to 1 by monte carlo integration 
println("Test LinProc.sample_p, LinProc.sample_p0, LinProc.lp, LinProc.lp0")

# monte carlo integration of the transition density p with proposals according to p0
function test1(h, x0, N, B, beta, A, lambda)
	su = 0
	#mu = x0
	mu =  B*x0+beta
	gamma = inv(A)
	l = chol(A)
	for n in [1:N]
		y = LinProc.sample_p0(h, x0,mu, l )
		su +=  exp(LinProc.lp(h, x0, y, B, beta, lambda)-LinProc.lp0(h, x0, y, mu, gamma))
	end
	su /= N
end	


# monte carlo integration of the transition density p0 with proposals according to p
function test2(h, x0, N, B, beta, A, lambda)
	su = 0
	mu =  B*x0+beta
	gamma = inv(A)

	for n in [1:N]
		y = LinProc.sample_p(h, x0, B, beta, lambda) 
		su +=  exp(LinProc.lp0(h, x0, y, mu, gamma)-LinProc.lp(h, x0, y, B, beta, lambda))
	end
	su /= N
end	

#wrapper for both tests

function test0(N)
	x0  = [0.2,0.2]
	A = [[1 0.2],[0.2 1]]
	B = [[-0.1 0.3],[ -0.3 -0.1]]
	h = 0.2
	beta = [0.2, 0.1]
	mu =  B*x0+beta

	lambda = Lyap.lyap(B', -A)
	gamma = inv(A)
	l = chol(A)
	su1 = test1(h, x0, N, B, beta, A, lambda )
#	println("Integral p(t,x; T, y, B, beta, a):", su1)

	su2 = test2(h, x0, N, B, beta, A, lambda )
#	println("Integral p0(t,x; T, y, B, beta, a):", su2)
	[su1, su2]
end
srand(5)
n1, n2 = test0(1E3)
#println("Check transition density p, p0 and sampler for p0, p")
#println("using seeded monte carlo integration.")
@test norm(n1-1.) < 5E-3
@test norm(n2-1.) < 5E-3
