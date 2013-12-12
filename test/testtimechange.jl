# testing clark.jl, takes some while and gives qualitative results, so not included in tests


srand(3)
d = 2

T = 1.7
t = 0.3

Broh = [ -2.21225   -0.49590;   0.631753  -1.34795]
A = .2*[ 0.170252  0.178533; 0.178533  2.43255]
 	

beta = [0.5,-0.5]
u = [1.,0.]
v = [0.5, -0.2]

function b(t,x) exp(-0.2*t)*Broh*x + beta end
B = 1.1exp(-0.2*T)*Broh

function ph(t,s) 5.*exp(-0.2*s)-5.*exp(-0.2*t) end
SIG = sqrtm(A)
sigma(t,x) = exp(-0.1*(T-t))*SIG
a(t,x) = exp(-0.2*(T-t))*A

Phi(t,s) = expm(ph(t,s)*Broh)


lambda = Lyap.lyap(B', -a(T,v))

 
N = 50

x = rand(d)
s = LinProc.tau(t,T)

 



Ts = linspace(0., T, N)
Dt = diff(Ts)
dt = Dt[1]
Smax = T
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
	ll =  LinProc.llikeliXcirc(0., T, Y, b, a, B, beta, lambda)
	Yt0[:,i] = Y[:,nt]
	Yll[i] = ll
end	

println("X° at ", round(t0,3), "($nt): ", "[", mean(Yt0[1,:]), ", ",mean(Yt0[2,:]), "]"	)


x1 = zeros(2, M)
s0 = LinProc.taui(t0, T)
ns = int(floor(s0/ds))
s0 = ds*ns

u0 =  LinProc.UofX(0,u,  T, v,  B, beta)
for i in 1:M
	DW = randn(2, N-1) .* sqrt(ds)
	U = LinProc.eulerv(0.0, u0, LinProc.bU(T, v, b, a, B, beta, lambda), (s,x) -> -sqrt(2/(T*(T-s))) *sigma(LinProc.ddd(s,x, 0., T, v,  B, beta)...), Ds, DW)
 	ll = LinProc.llikeliU(S, U, 0., T, v, b, a,  B, beta, lambda)
	Y2t0[:,i] =  LinProc.XofU(s0, U[:,ns],  T, v,  B, beta)
	Y2ll[i] = ll
end

println("X°(U) at ", round(LinProc.tau(s0,T),3), "    : ", "[",mean(Y2t0[1,:]), ",", mean(Y2t0[2,:]), "]")
	 			 

println("T-dt = ", round(T-dt,3), " ~ ", round(LinProc.tau(T-ds,T),3))




plambda = exp(LinProc.lp(T, u, v, B, beta, lambda))
p = exp(lp(0., T, u,v, ph, Broh, beta, s -> a(s,NaN)))
#println(mc(exp(imag(Xll))), mc(exp(imag(X2ll)))) 

println("~p ", round(plambda,5), " p ", round(p,5), " ", repr(mc(exp(real(Yll))*plambda)), repr(mc(exp(real(Y2ll))*plambda)) )
w1 = mc(plambda/p*exp(real(Yll)))
w2 = mc(plambda/p*exp(real(Y2ll)))
println("Check Σ ptilda/p*psi(X°) = ", repr(w1), " ≈ 1", " and Σ ptilda/p*psi(U)", repr(w2), " ≈ 1")
#Test.@test abs(w1[1] - 1.) < 1.96w1[2]
#Test.@test abs(w2[1] - 1.) < 1.96w2[2]


