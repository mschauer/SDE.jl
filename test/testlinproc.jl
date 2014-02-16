using SDE
using Base.Test
require(Pkg.dir("SDE","src", "Randm.jl"))
#include(Pkg.dir("SDE","src", "quad.jl"))
import SDE.eps2

for d in [1 2 3]
    xi = 0.3
	global T, t, v, h, P, B, A, beta, s
	println("Test dimension $d xi $xi")

	if (d== 1) srand(7)
	elseif (d== 2) srand(8)
	else srand(5)
	end
	B = xi.*Randm.randstable(d)
	A = Randm.randposdef(d)
	Si = chol(A, :L)
	beta = 0.5*randn(d)
	
	P = MvLinPro(B, beta, Si)
	
	function aa(s)
	 A
	end
	

	T = 1.7
	t = 0.5
	v = zeros(d*d)



	println("Test SDE.H, SDE.K")

	H1 = SDE.H(t,T, P)
	H2 = SDE.varH(t, T, B, aa)
	K1 = SDE.K(t, T, P)
	H3 = expm((T-t)*B')*inv(K1)*expm((T-t)*B)


	h = T-t
	phi = expm(h*B)
	phim = expm(-h*B)
	K2 = similar(H1)
	H4 = similar(H1)
	K2 = P.lambda - phi*P.lambda*phi'
	H4 = phi'*inv(P.lambda - phi*P.lambda*phi')*phi
	
	@test norm(H1 - H2) < 1E-10
	@test norm(H1 - H3) < 1E-10
	@test norm(H1 - H4) < 1E-10

	@test norm(K1-K2) < d*24*eps()

	function varr3(t, x, T, v, P)
	    h = T - t
		binv = inv(P.B)
		phi = expm(h*P.B)
		phim = expm(-h*P.B)
		vt =  phim*( v) + phim * binv * beta -binv*beta 
		SDE.H(t, T, P) * (vt-x )
	
	end

	println("Test SDE.r")

	srand(5)
	x0 = x = randn(d)
	v = randn(d)
	r1 = SDE.r(t, x, T, v, P)
	r2 = SDE.varr(t, x, T, v, P)
	@test norm(r1 - r2) < 1E-10

	if(norm(B) > eps2)
		r3 = varr3(t, x, T, v, P)
		@test norm(r1 - r3) < 1E-10
	end

	println("Test SDE.mu")

	mu1 = SDE.mu(t, x, T, P)
	mu2 = SDE.varmu(t, x, T, P)
	@test norm(mu1 - mu2) < 1E-13


	#println(mu1)
	#println(mu2)

	#checks if transition density p integrates to 1 by monte carlo integration 
	println("Test SDE.sample_p, SDE.sample_p0, SDE.lp, SDE.lp0")

	# monte carlo integration of the transition density p with proposals according to p0
	# monte carlo integration of the transition density p0 with proposals according to p
	function testp(h, x0, N, B, beta, Si)
		mu =  B*x0+beta
		P = MvLinPro(B, beta, Si)
		P0 = MvAffPro(mu, Si)
		su1 = 0.
		for n in 1:N
			y = SDE.samplep(0., x0, h, P0 )
			su1 +=  exp(SDE.lp(0., x0, h,  y, P)-SDE.lp(0., x0, h, y, P0))
		end
		su1 /= N
		
		su2 = 0.
		for n in 1:N
			y = SDE.samplep(0., x0, h, P )
			su2 +=  exp(SDE.lp(0., x0, h,  y, P0)-SDE.lp(0., x0, h, y, P))
		end
		su2 /= N
		return su1, su2
	end	


	srand(5)
	n1, n2 = testp(0.2/d, x0, 1E4, B, beta, Si)
	#println("Check transition density p, p0 and sampler for p0, p")
	#println("using seeded monte carlo integration.")
	println(n1, " ", n2)
	@test norm(n1-1.) < 5E-2
	@test norm(n2-1.) < 5E-2

	@test norm(SDE.lp(t, x,T, v, P) - SDE.varlp(t, x,T, v, (t,s) -> t-s, B, beta, A)) < 1E-10



	s = SDE.soft(t, 0., T)

	v1 = SDE.Vs(s, T, v, P)
	v2 = SDE.V(t, T, v, P)
	
	@test norm(v1-v2) < d^2*1E-13

	us =  SDE.uofx(s, x, T, v,  P)
	xt =  SDE.xofu(s, us,  T, v,  P)
	@test norm(x - xt) < 1E-10

	r1 =  SDE.r(t, x, T, v, P)
 	r2 =  SDE.J(s, T, P) *us*T/(T-s)
	@test norm(r1 - r2) < 1E-10

end

