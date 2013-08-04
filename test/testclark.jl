# test clark.jl

include("../src/linproc.jl")
include("../src/randm.jl")
include("../src/lyap.jl")
include("../src/quad.jl")
require("misc.jl")

using Base.Test
srand(5)
d = 2

T = 1.7
t = 0.5

B = [ -2.21225   -0.49590;   0.631753  -1.34795]
A = 5*[ 0.170252  0.178533; 0.178533  2.43255]
lambdal = Lyap.lyap(B', -A)


beta = [0.5,-0.5]
u = [1.,0.]

b(t,x) = exp(-0.2*t)*B*x + beta
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
SIG = sqrtm(A)
sigma(t,x) = exp(-0.1*(T-t))*SIG
a(t,x) = exp(-0.2*(T-t))*A

 
N = 100



v = [0.5, -0.1]
lambda = Lyap.lyap(B', -a(T,v))

Ts = linspace(0., T, N)
Dt = diff(Ts)
dt = Dt[1]
Smax = LinProc.taui(T-dt,T)
S = linspace(0.,Smax, N)
Ds = diff(S)
ds = Ds[1]

M = 2000 

Yll = zeros(M) 
Y2ll = zeros(M) 
Y3ll = zeros(M) 

Yt0 = zeros(2, M)
Y2t0 = zeros(2, M)
 
nt = int(floor(5N/10))

t0 =  dt*nt

for i in 1:M
	DW = randn(2, N-1) .* sqrt(dt)
	Y = LinProc.eulerv(0.0, u, v, LinProc.bcirc(T, v, b, sigma, B, beta, lambda), sigma, Dt, DW)
	ll =  LinProc.llikeliXcirc(0, T, Y, b, a, B, beta, lambda)
	Yt0[:,i] = Y[:,nt]
	Yll[i] = ll
end	

println("X° at ", round(t0,3), "($nt): ", [mean(Yt0[1,:]), mean(Yt0[2,:])])


x1 = zeros(2, M)
s0 = LinProc.taui(t0, T)
ns = int(floor(s0/ds))
s0 = ds*ns
 
for i in 1:M
	DW = randn(2, N-1) .* sqrt(ds)
	u0 =  LinProc.UofX(0,u,  T, v,  B, beta)
	U = LinProc.eulerv(0.0, u0, LinProc.bU(T, v, b, a, B, beta, lambda), (s,x) -> sqrt(T)*sigma(LinProc.ddd(s,x, 0., T, v,  B, beta)...), Ds, DW)
 	ll = LinProc.llikeliU(S, U, 0., T, v, b, a,  B, beta, lambda)
	Y2t0[:,i] =  LinProc.XofU(s0, U[:,ns],  T, v,  B, beta)
	Y2ll[i] = ll
end

println("X°(U) at ", round(LinProc.tau(s0,T),3), "    : ", [mean(Y2t0[1,:]), mean(Y2t0[2,:])])

S = -log(1.-Ts/T) 
Ds = diff(S)
	 	 
 
		
for i in 1:M
	DW = randn(2,N-1).*[sqrt(Ds)  sqrt(Ds)]'
	u0 =  LinProc.UofX(0,u,  T, v,  B, beta)
	U = LinProc.eulerv(0.0, u0, LinProc.bU(T, v, b, a, B, beta, lambda), (s,x) -> sqrt(T)*sigma(LinProc.ddd(s,x, 0., T, v,  B, beta)...), Ds, DW)
 	ll = LinProc.llikeliU(S, U, 0., T, v, b, a,  B, beta, lambda)
	Y3ll[i] = ll
end


	 			 

println("T-dt = ", round(T-dt,3), " ~ ", round(LinProc.tau(N*ds,T),3))

println(mc(exp(real(Yll))), mc(exp(real(Y2ll)))) 
#println(mc(exp(imag(Xll))), mc(exp(imag(X2ll)))) 

plambda = exp(LinProc.lp(T, u, v, B, beta, lambda))
p = exp(lp(0, T, u,v, ph, B, beta, s -> a(s,NaN)))

println("~p ", round(plambda,5), " p ", round(p,5), " ", mc(exp(real(Yll))*plambda), mc(exp(real(Y2ll))*plambda) )
println()


