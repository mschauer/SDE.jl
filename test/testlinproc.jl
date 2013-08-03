# test linproc.jl

include("../src/linproc.jl")
include("../src/randm.jl")
include("../src/lyap.jl")

using Cubature
using Base.Test


d= 2
# B = [-10.3268 .67701   -1.85162;
#   2.60125  -5.35212    0.4388; 
#  -1.42184   0.904901  -0.5423]
#
srand(5)
B = Randm.randstable(d)
A = Randm.randposdef(d)
function aa(s)
 A
end
lambdal = Lyap.lyap(B', -A)
lambdas = Lyap.syl(B, B',-A)
lambda = lambdal

T = 1.7
t = 0.5
v = zeros(d*d)

# using Cubature.jl
function Hquad(t, T, B, a::Function)
	function f(s, v)
		v[:] = vec(expm(-(s-t)*B)*a(s)*expm(-(s-t)*B)')
	end
	Q = reshape(hquadrature(d*d, f, t, T, 1E-15, 1E-15, 5000)[1], d,d)
	inv(Q)
end	


println("Test LinProc.H, LinProc.K")

H1 = LinProc.H(T-t, B, lambda)
H2 = Hquad(t, T, B, aa)
K1 = LinProc.K(T-t, B, lambda)
H3 = expm((T-t)*B')*inv(K1)*expm((T-t)*B)


h = T-t
phi = expm(h*B)
phim = expm(-h*B)
K2 = lambda - phi*lambda*phi'
H4 = phi'*inv(lambda - phi*lambda*phi')*phi
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
	LinProc.H(h, B, lambda) * (vt-x )
	
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

s = -log(1-t/T) 
v1 = LinProc.Vtau (s,T, v, B, beta)
v2 = LinProc.V(T-t, v, B, beta)

@test norm(v1-v2) < 1E-13
v3 = LinProc.Vtau (27,T, v, B, beta)
@test norm(v-v3) < 1E-10

J1 = T*exp(-s)*LinProc.H(T-t, B, lambda)
J2 = LinProc.J1(s,T, B, A, lambda)
@test norm(J1-J2) < 1E-10
#@test isposdef(J2)

us =  LinProc.UofX(s,x,  T, v,  B, beta)
xt =  LinProc.XofU(s,us,  T, v,  B, beta)
@test norm(x - xt) < 1E-10


if (d == 2)
	a1 = LinProc.J1(6,T, B, A, lambda)
	a2 = LinProc.J2(6,T, B, A, lambda)
	@test norm(a1-a2) < 1E-10
	a3 = LinProc.J2(27,T, B,A, lambda)
 	@test norm(a3- inv(A)) < 1E-10

end


#julia> A
#2x2 Float64 Array:
# 0.170252  0.178533
# 0.178533  2.43255

#julia> B
#2x2 Float64 Array:
# -2.21225   -0.49590
#  0.631753  -1.34795

beta = [0.5,-0.5]
u = [1.,0.]
A = 5*[ 0.170252  0.178533; 0.178533  2.43255]

#b(t,x) = exp(-0.2*t)*B*x + beta
b(t,x) = exp(-0.2*t)*B*x + beta
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
SIG = sqrtm(A)
sigma(t,x) = exp(-0.1*(T-t))*SIG
a(t,x) = exp(-0.2*(T-t))*A

 
#ph(t,s) = t-s
function Qquad(s, T, ph, B, a::Function)
	
	function f(tau, v)
		v[:] = vec(expm(ph(s,tau)*B)*a(tau)*expm(ph(s,tau)*B)')
	end
	Q = reshape(hquadrature(d*d, f, s, T, 1E-15, 1E-15, 5000)[1], d,d)
	Q
end
 
function Vquad(s,T, v, ph, b, beta)
	function f(tau, y)
		y[:] = expm(ph(s,tau)*b)*beta(tau)
	end
	expm(-ph(T,s)*b)*v - hquadrature(d, f, s, T, 1E-15, 1E-15, 5000)[1]
	
end


function lp(t, T, x, y, ph, B, beta, a::Function)
 	z = (x -  Vquad(t,T, y, ph, B, t -> beta))
	Q = Qquad(t, T, ph,  B, a)
	l = chol(Q, :L)
	K =  expm(ph(T,t)*B)*Q*expm(ph(T,t)*B)'
	(-1/2*length(x)*log(2pi) -log(apply(*,diag(chol(K)))) - 0.5*norm(l\z)^2) #  - 0.5*log(det(K(h,b, lambda)))
end
#println("LinProc.lp ", LinProc.lp(T-t, u, v, B, beta,  Lyap.lyap(B', -A)), " lp ", lp(t, T, u,v, (t,s) -> t-s, B, beta, s -> A))

@test norm(LinProc.lp(T-t, u, v, B, beta,  Lyap.lyap(B', -A)) - lp(t, T, u,v, (t,s) -> t-s, B, beta, s -> A)) < 1E-10

N = 50


require("misc.jl")
Dt = diff(linspace(0., T, N))
dt = Dt[1]
#DW = randn(2, N-1) .* sqrt(dt)
#yy = LinProc.eulerv(0.0, u,  b(B,beta), sigma(A), Dt, DW)
#v = yy[:,N]


v = [0.5, -0.1]
lambda = Lyap.lyap(B', -a(T,v))


Dt = diff(linspace(0., T, N))
dt = Dt[1]
Smax = LinProc.taui(T-dt,T)
S = linspace(0.,Smax, N)
Ds = diff(S)
ds = Ds[1]
yy = zeros(2, N)

M = 20000
Xt0 = zeros(2, M)
XT = zeros(2, M)
Xll = zeros(M) 
X2ll = zeros(M) 
X3ll = zeros(M) 

X2t0 = zeros(2, M)
Ut0 = zeros(2, M)
nt = floor(5N/10)

t0 =  dt*nt

#for i in 1:M
#	DW = randn(2, N-1) .* sqrt(dt)
#	yy = LinProc.eulerv(0.0, u, b, sigma, Dt, DW)
# 	XT[:,i] = yy[:,end]
#	 
#end
#println(mean(XT[1,:]), mean(XT[2,:]))

 

for i in 1:M
	DW = randn(2, N-1) .* sqrt(dt)
	yy = LinProc.eulerv(0.0, u, v, LinProc.Bcirc(T, v, b, sigma, B, beta, lambda), sigma, Dt, DW)
	ll =  LinProc.llikelixcirc(0, T, yy, b, a, B, beta, lambda)
	Xt0[:,i] = yy[:,nt]
	Xll[i] = ll
end	

println("X° at ", round(t0,3), "($nt): ", mean(Xt0[1,:]), " ", mean(Xt0[2,:]))


x1 = zeros(2, M)
s0 = LinProc.taui(t0, T)
ns = floor(s0/ds)
s0 = ds*ns
 
for i in 1:M
	DW = randn(2, N-1) .* sqrt(ds)
	u0 =  LinProc.UofX(0,u,  T, v,  B, beta)
	yy = LinProc.eulerv(0.0, u0, LinProc.bU(T, v, b, a, B, beta, lambda), (s,x) -> sqrt(T)*sigma(LinProc.ddd(s,x, 0., T, v,  B, beta)...), Ds, DW)
#	println(size(yy))
	ll = LinProc.llikeliU(S, yy, T, v, b, a,  B, beta, lambda)
	Ut0[:,i] = yy[:,ns]
	X2t0[:,i] =  LinProc.XofU(s0, yy[:,ns],  T, v,  B, beta)
	X2ll[i] = ll
end


S = -log(1.-Ts/T) 
Ds = diff(S)
	 	 
 
		
for i in 1:M
	DW = randn(2,n-1).*[sqrt(Ds)  sqrt(Ds)]'
	u0 =  LinProc.UofX(0,u,  T, v,  B, beta)
	yy = LinProc.eulerv(0.0, u0, LinProc.bU(T, v, b, a, B, beta, lambda), (s,x) -> sqrt(T)*sigma(LinProc.ddd(s,x, 0., T, v,  B, beta)...), Ds, DW)
#	println(size(yy))
	ll = LinProc.llikeliU(S, yy, T, v, b, a,  B, beta, lambda)
	X3ll[i] = ll
end


	 			
println("U at ", round(s0,3), "($ns): ", mean(Ut0[1,:]), " ", mean(Ut0[2,:]))
println("X°° at ", round(LinProc.tau(s0,T),3), "    : ", mean(X2t0[1,:]), " ", mean(X2t0[2,:]))

println("T-dt = ", round(T-dt,3), " ~ ", round(LinProc.tau(N*ds,T),3))

println(mc(exp(real(Xll))), mc(exp(real(X2ll)))) 
#println(mc(exp(imag(Xll))), mc(exp(imag(X2ll)))) 

plambda = exp(LinProc.lp(T, u, v, B, beta, lambda))
p = exp(lp(0, T, u,v, ph, B, beta, s -> a(s,NaN)))

println("~p ", round(plambda,5), " p ", round(p,5), " ", mc(exp(real(Xll))*plambda), mc(exp(real(X2ll))*plambda) )
println()


