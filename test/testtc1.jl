using SDE
using Cubature
import SDE.b, SDE.sigma, SDE.a

#srand(7)
if !isdefined(:Model)

type Model <: UvPro
    B
    beta
    Sigma
    A
    d::Int   
end
end

b(t,x, P::Model) =   exp(-0.2*t)*P.B*x + P.beta
sigma(t,x, P::Model) = exp(-0.1*(t))*P.Sigma
a(t,x,  P::Model) = exp(-0.2*(t))*P.A

let d = 1

Q95 = sqrt(2.)*erfinv(0.95)
function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, Q95*ste], int(2-log(10,Q95*ste)))
end

T = tmax = 1.7
tmin = 0.
 
B0 = -2.21225
beta0 = 0.1
sigma0 = .15
ph0(t,s ) =  5.*exp(-0.2*s)-5.*exp(-0.2*t)

Model() = Model(B0, beta0,sigma0, sigma0^2, 1)


u = 0.3
v = 0.5

#taken slightly off 
B = exp(-0.2*T)*B0 
 
N = 200
 

P = Model()
Pt = UvLinPro(B, beta0, sigma(T, v, P))


tt = linspace(0., T, N)
dt = (tt[N]-tt[1])/(N-1)

ss = linspace(0.,T, N)
ds = (ss[N]-ss[1])/(N-1)

ttofss = map(s -> tofs(s, 0., T), ss)

M = 20000

Yll = zeros(M) 
UYll = zeros(M)

Y2ll = zeros(M) 
UY2ll = zeros(M) 

Ull = zeros(M) 
YUll = zeros(M)

Yt0 = zeros(M)
Y2t0 = zeros(M)
YUt0 = zeros(M)
Us0 = zeros(M)
 
nt = int(floor(5N/10))


t0 =  dt*(nt-1)
W = sample(ss, Wiener())

@time for i in 1:M
	resample!(W, Wiener())
	Y = guidedeuler(u, W, T, v, Pt, P)
	U = UofX(Y, 0., T, v, Pt) 
	ll = llikeliXcirc(Y, Pt, P)
	ll2 = llikeliU(U, 0., T, v, Pt, P)
	Yt0[i] = Y.yy[nt]
	Yll[i] = ll
	UYll[i] = ll2
end	

println("X° at ", round(t0,3), "($nt): ""[", repr(mc(Yt0)),"]")



x1 = zeros(M)
s0 = soft(t0, 0, T)
ns = 1 + int(floor(s0/ds))
s0 = ds*(ns-1)



u0 = uofx(0.0, u,  T, v, Pt)
W = sample(ss, Wiener())
U = eulerU(u0, W, 0., T, v, Pt, P)
@time for i in 1:M
 	resample!(W, Wiener())
	eulerU!(U, u0, W, 0., T, v, Pt, P)
	YU = XofU(U, 0., T, v, Pt) 
	llu = llikeliU(U, 0., T, v, Pt, P)
	llyu = llikeliXcirc(YU, Pt, P)
	Us0[i] = U.yy[ns]
	YUt0[i] = xofu(s0, U.yy[ns],  T, v, Pt)
	Ull[i] = llu
	YUll[i] = llyu
end



println("X°(U) at ", round(tofs(s0,0., T),3), "($ns): ",  "[", repr(mc(YUt0[:])),"]")

 
W = sample(ttofss, Wiener())
@time for i in 1:M
	resample!(W, Wiener())
	Y2 = guidedeuler(u, W, T, v, Pt, P)
	Y2ll[i]  = llikeliXcirc(Y2, Pt, P)
	Y2t0[i] = Y2.yy[ns]
	U2 = UofX(Y2, 0., T, v, Pt) 
	UY2ll[i] = llikeliU(U2, 0., T, v, Pt, P)
    
end	

println("X°° at ", round(ttofss[ns],3), "($ns): ", "[", repr(mc(Y2t0[:])),  "]"	)
	 			 



    pt = exp(lp(0., u, T, v, Pt))
	p = exp(SDE.varlp(0., u, T, v, ph0, B0, beta0, s->a(s, u, P)))
	println("p ", round(p,5), " ", repr(mc(exp(Yll)*pt)), repr(mc(exp(Ull)*pt)) , repr(mc(exp(Y2ll)*pt)))
	println("~p ", round(pt,5))
	w1 = mc(pt/p*exp(Yll))
	w1b = mc(pt/p*exp(UYll))
	
	w2 = mc(pt/p*exp(Ull))
	w2b = mc(pt/p*exp(YUll))
	w3 = mc(pt/p*exp(Y2ll))
    w3b = mc(pt/p*exp(UY2ll))	
    
	println(" Σ ~p/p*psi(X°) = ", repr(w1), " ≈ 1\n",
	        " Σ ~p/p*psi(UX°) = ", repr(w1b), " ≈ 1\n",
            " Σ ~p/p*psi(U) = ", repr(w2), " ≈ 1\n",
            " Σ ~p/p*psi(YU) = ", repr(w2b), " ≈ 1\n", 
            " Σ ~p/p*psi(Xs)= ", repr(w3), " ≈ 1\n",
            " Σ ~p/p*psi(UXs)= ", repr(w3b), " ≈ 1"
            )
	println("lmax ", maximum(exp(Yll)), " ", maximum(exp(Ull)), " ", maximum(exp(Y2ll)))

    Test.@test abs(w1[1] - 1.) < 1.96w1[2]
    Test.@test abs(w2[1] - 1.) < 1.96w2[2]
    Test.@test abs(w3[1] - 1.) < 1.96w3[2]

end

