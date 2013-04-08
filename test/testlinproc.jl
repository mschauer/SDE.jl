# test linproc.jl

using Cubature
using LinProc
using Randm
using Lyap
using Test

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
println(v1)

mu1 = mu(h, x, B, beta)
mu2 = varmu(h, x, B, beta)
# inv(phim*(lambda - phi*lambda*phi')*phim')^C

#println(mu1)
#println(mu2)
