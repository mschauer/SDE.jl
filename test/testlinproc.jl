# test linproc.jl
require("SDE")

using Cubature
using LinProc
using Randm
using Lyap
using Base.Test

d= 4
# B = [-10.3268 .67701   -1.85162;
#   2.60125  -5.35212    0.4388; 
#  -1.42184   0.904901  -0.5423]
#
B = randstable(d)
A = randposdef(d)
function a(s)
 A
end
lambda = lyap(B', -A)

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
H1 = LinProc.H(T-t, B, lambda)
H2 = Hquad(t, T, B, a)
K1 = LinProc.K(T-t, B, lambda)
H3 = -expm((T-t)*B')*inv(K1)*expm((T-t)*B)


h = T-t
phi = expm(h*B)
phim = expm(-h*B)
K2 = lambda - phi*lambda*phi'
H4 = -phi'*inv(lambda - phi*lambda*phi')*phi
H5 = LinProc.H(h, B, lambda)

@test norm(H1 - H2) < 1E-5
@test norm(H3 - H2) < 1E-5

#obtaining r via quadrature
function varmu(h, x, B, beta)
	function f(s, y)
		y[:] = expm(-s*B)*beta
	end
	integral = hquadrature(d, f, 0, h, 1E-15, 1E-15, 800)[1]
	expm(h*B)*(x + integral)
end	

function mu(h, x, B, beta)
	
	binv = inv(B)
	phi = expm(h*B)
	phim = expm(-h*B)
	integral = binv*beta - phim * binv * beta
	phi*(x + integral)
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


x = randn(d)
v = randn(d)
beta = randn(d)
r1 = LinProc.r(T-t, x, v, B, beta, lambda)
r2 = varr(T-t, x, v, B, beta, lambda)
r3 = varr3(T-t, x, v, B, beta, lambda)

println(r1)
println(r2)
println(r3)
println(v)

mu1 = mu(h, x, B, beta)
mu2 = varmu(h, x, B, beta)
# inv(phim*(lambda - phi*lambda*phi')*phim')^C

#println(mu1)
#println(mu2)

#checks if transition density p integrates to 1 by monte carlo integration 
function test1(h, x0, N, B, beta, A, lambda)
	su = 0
	mu =  B*x0+beta
	mu = x0
	gamma = inv(A)
	l = chol(A)
	for n in [1:N]
		y = sample_p0(h, x0,mu, l )
		#su +=  exp(lp0(h, x0, y, mu, diagm([1,1])) +   0*lp(h, x0, y, B, beta, lambda)-lp0(h, x0, y, mu, gamma))
		su +=  exp(lp(h, x0, y, B, beta, lambda)-lp0(h, x0, y, mu, gamma))
	end
	su /= N
end	

function test0()
 x0  = [0.2,0.2]
 A = [[1 0.2],[0.2 1]]
 B = [[-0.1 0.3],[ -0.3 -0.1]]
 h = 0.2
 beta = [0.2, 0.1]
 mu =  B*x0+beta

 lambda = lyap(B', -A)
 gamma = inv(A)
 l = chol(A)
 su = LinProc.test1(h, x0, 1E5, B, beta, A, lambda )

 println("Integral p(t,x; T, y, B, beta, a):", su)


end

