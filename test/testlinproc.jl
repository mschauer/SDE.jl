include(Pkg.dir("SDE","src", "LinProc.jl"))
require(Pkg.dir("SDE","src", "Randm.jl"))


for d in [1 2 3]
for xi in [0.0 0.5]
	global T, t, v, h, lambda, B, A, beta, s
	println("Test dimension $d xi $xi")
	# B = [-10.3268 .67701   -1.85162;
	#   2.60125  -5.35212    0.4388; 
	#  -1.42184   0.904901  -0.5423]
	#
	if (d== 1) srand(7)
	elseif (d== 2) srand(8)
	else srand(5)
	end
	B = xi.*Randm.randstable(d)
	A = Randm.randposdef(d)
	function aa(s)
	 A
	end
	if norm(B) > eps()
		lambdal = SDE.lyap(B', -A)
		lambdas = SDE.syl(B, B',-A)
	else
		lambdal = inv(A)
		lambdas = inv(A)
	end

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
	K2 = similar(H1)
	H4 = similar(H1)
	if (norm(B) < eps2)
		K2 = A*h
		H4 = phi'*lambda/h*phi

	else 
		K2 = lambda - phi*lambda*phi'
		H4 = phi'*inv(lambda - phi*lambda*phi')*phi

	end
	
	H5 = LinProc.H(h, B, lambda)

	@test norm(H1 - H2) < 1E-10
	@test norm(H1 - H3) < 1E-10
	@test norm(H1 - H4) < 1E-10
	@test norm(H1 - H5) < 1E-10

	@test norm(K1-K2) < d*24*eps()

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
	beta = 0.5*randn(d)
	r1 = LinProc.r(T-t, x, v, B, beta, lambda)
	r2 = varr(T-t, x, v, B, beta, lambda)
	@test norm(r1 - r2) < 1E-10

	if(norm(B) > eps2)
		r3 = varr3(T-t, x, v, B, beta, lambda)
		@test norm(r1 - r3) < 1E-10
	end

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
		h = 0.1
		mu =  B*x0+beta
		lambda = gamma = inv(A)
		if (norm(B) > eps2)
			lambda = SDE.lyap(B', -A)
		end
		
		l = chol(A)
		su1 = test1(h/2, x0, N, B, beta, A, lambda )
	#	println("Integral p(t,x; T, y, B, beta, a):", su1)

		su2 = test2(h/2, x0, N, B, beta, A, lambda )
	#	println("Integral p0(t,x; T, y, B, beta, a):", su2)
		[su1, su2]
	end
	srand(5)
	n1, n2 = test0(1E4, B)
	#println("Check transition density p, p0 and sampler for p0, p")
	#println("using seeded monte carlo integration.")
	println(n1, " ", n2)
	@test norm(n1-1.) < 5E-2
	@test norm(n2-1.) < 5E-2

	@test norm(LinProc.lp(T-t, x, v, B, beta,  lambda) - lp(t, T, x,v, (t,s) -> t-s, B, beta, s -> A)) < 1E-10


#	b(t,x) = exp(-0.2*t)*B*x + beta
#	ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
#	SIG = sqrtm(A)
#	sigma(t,x) = exp(-0.1*(T-t))*SIG
#	a(t,x) = exp(-0.2*(T-t))*A


	s = LinProc.taui(t, T)

	v1 = LinProc.Vtau (s,T, v, B, beta)
	v2 = LinProc.V(T-t, v, B, beta)
	
	@test norm(v1-v2) < d^2*1E-13
	#for clark
#	v3 = LinProc.Clark.Vtau (27,T, v, B, beta) 
#	@test norm(v-v3) < 1E-10


	v3 = LinProc.Vtau (T,T, v, B, beta)
	@test norm(v-v3) < 1E-10

#	J1 = T*exp(-s)*LinProc.H(T-t, B, lambda)
	J1 = (T-s)^2/T*LinProc.H(T-t, B, lambda)
	J2 = LinProc.J(s,T, B, A, lambda)
	println(norm(J1-J2))

	@test norm(J1-J2) < 1E-10
	us =  LinProc.UofX(s,x,  T, v,  B, beta)
	xt =  LinProc.XofU(s,us,  T, v,  B, beta)
	@test norm(x - xt) < 1E-10


	r1 =  LinProc.H(T-LinProc.tau(t,T), B, lambda)*(LinProc.V(T-LinProc.tau(t,T), v, B, beta)-LinProc.XofU(t,x,  T, v,  B, beta))
	r2 =  LinProc.H(T-LinProc.tau(t,T), B, lambda)*(T-t)*x
 	r3 =  LinProc.J(t, T, B, A, lambda) *x*T/(T-t)
	r4 = LinProc.r(T-LinProc.tau(t,T), LinProc.XofU(t,x,  T, v,  B, beta), v, B, beta, lambda)
	@test norm(r3 - r4) < 1E-10

	print(LinProc.J(t, T, B, A, lambda) *x*T/(T-t))
	

end
end
