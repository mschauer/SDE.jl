# test linproc.jl

include("../src/linproc.jl")
include("../src/randm.jl")
include("../src/lyap.jl")
include("../src/quad.jl")

using Base.Test

for d in [1 2 3]
	println("Test dimension =$d")
	# B = [-10.3268 .67701   -1.85162;
	#   2.60125  -5.35212    0.4388; 
	#  -1.42184   0.904901  -0.5423]
	#
	if (d== 1) srand(7)
	elseif (d== 2) srand(8)
	else srand(5)
	end
	B = 0.5*Randm.randstable(d)
	A = Randm.randposdef(d)
	function aa(s)
	 A
	end
	lambdal = Lyap.lyap(B', -A)
	lambdas = Lyap.syl(B, B',-A)

	@test norm(lambdal - lambdas) < 1E-10
	lambda = lambdal

	T = 1.7
	t = 0.5
	v = zeros(d*d)



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
	x0 = x = randn(d)
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
	@test norm(mu1 - mu2) < 1E-13


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

	function test0(N, B)
	#	x0  = [0.2,0.2]
	#	A = [[1 0.2],[0.2 1]]
	#	B = [[-0.1 0.3],[ -0.3 -0.1]]
		h = 0.1
	#	beta = [0.2, 0.1]
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
	n1, n2 = test0(5E3, 0.5B)
	#println("Check transition density p, p0 and sampler for p0, p")
	#println("using seeded monte carlo integration.")
	println(n1, " ", n2)
	@test norm(n1-1.) < 5E-2
	@test norm(n2-1.) < 5E-2

	u = x

	b(t,x) = exp(-0.2*t)*B*x + beta
	ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
	SIG = sqrtm(A)
	sigma(t,x) = exp(-0.1*(T-t))*SIG
	a(t,x) = exp(-0.2*(T-t))*A

	@test norm(LinProc.lp(T-t, u, v, B, beta,  Lyap.lyap(B', -A)) - lp(t, T, u,v, (t,s) -> t-s, B, beta, s -> A)) < 1E-10

end
