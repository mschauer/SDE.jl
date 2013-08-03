# test linproc.jl

include("../src/linproc.jl")
include("../src/randm.jl")
include("../src/lyap.jl")
include("../src/quad.jl")
require("misc.jl")

using Base.Test

for d in [1, 3, 2]

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

srand(5)
x = randn(d)
v = randn(d)
beta = randn(d)

s = -log(1-t/T) 
v1 = LinProc.Vtau (s,T, v, B, beta)
v2 = LinProc.V(T-t, v, B, beta)

@test norm(v1-v2) < 1E-13
v3 = LinProc.Vtau (27,T, v, B, beta)
@test norm(v-v3) < 1E-10

J1 = T*exp(-s)*LinProc.H(T-t, B, lambda)
J2 = LinProc.J1(s,T, B, A, lambda)
@test norm(J1-J2) < 1E-10

us =  LinProc.UofX(s,x,  T, v,  B, beta)
xt =  LinProc.XofU(s,us,  T, v,  B, beta)
@test norm(x - xt) < 1E-10

end

D = 2

T = 1.7
t = 0.5

B = [ -2.21225   -0.49590;   0.631753  -1.34795]
A = 5*[ 0.170252  0.178533; 0.178533  2.43255]
lambdal = Lyap.lyap(B', -A)

a1 = LinProc.J1(6,T, B, A, lambda)
a2 = LinProc.J2(6,T, B, A, lambda)
@test norm(a1-a2) < 1E-10
a3 = LinProc.J2(27,T, B,A, lambda)
@test norm(a3- inv(A)) < 1E-10


beta = [0.5,-0.5]
u = [1.,0.]

b(t,x) = exp(-0.2*t)*B*x + beta
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
SIG = sqrtm(A)
sigma(t,x) = exp(-0.1*(T-t))*SIG
a(t,x) = exp(-0.2*(T-t))*A

 
N = 50



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


