using Lyap

require("linproc.jl")
require("misc.jl")

srand(3)


si = 0.5
# si = 0.05
d = 2

B = [ -2.21225   -0.49590;   0.631753  -1.34795]
A = [ 0.170252  0.178533; 0.178533  2.43255]

beta = [0.5,-0.5]

function sigma0(s,x) 
	m = norm(x)
	if (m > eps())  
		rho = 1 + 5*2atan(m)/pi
		return( si*((1-2atan(m)/pi)eye(2)+ 2atan(m)/pi/m*[[x[2], -x[1]]  [rho*x[1], rho*x[2]]]))
	end
	return si*eye(2)
end

function sigma1(s,x)
 	m = norm(x-u)
	[1 1; -1  (1 + 8atan(2m)/pi)]
end


#sigma(t,x) = sigma0(0., [sin(t), cos(t)])
#sigma(t,x) = sigma1(t,x)
sigma(t,x) = sigma0(t,x)




#sigma(t,x) = sigma0(t,x)


a(s,x) = sigma(s,x)*sigma(s,x)'



b(t,x) = exp(-0.2*t)*B*x + beta
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
#SIG = sqrtm(A)
#SIG = chol(A,:L)
#sigma(t,x) = exp(-0.1*(T-t))*SIG
#a(t,x) = exp(-0.2*(T-t))*A
u = [1.,0.]
 
T = 1.2
v0 = [0.5, -0.1]
v0max = v0 

# if v=re() is distributed with log density le()
# and yy is a path with drift LinProc.Bcirc(T, v, b, sigma, B, beta, lambda) and diffusion sigma
# (lambda has to be computed for B(T, v) and a(T,v)!)
# then the expectation of llikelicirc(... yy, ...) is 1/pbar
# where pbar = exp(lp(..., B, beta, lambda) - le(v0)))
		
 
function rare(K, N, E, re, pe)

	println("Compute P(X in E), K=$K")
 
	llmax = [Inf,-Inf]
	lplambdamax = [Inf,-Inf]

	L = L2 = 0.0
	V = 0.0
	V2 = 0.0

	# avoid cancellation 
	V0 = BigFloat(0.0)
	V02 = BigFloat(0.0)
	
	Dt = diff(linspace(0., T, N))
	dt = Dt[1]
	Smax = LinProc.taui(T-dt,T)
	S = linspace(0.,Smax, N)
	Ds = diff(S)
	ds = Ds[1]

	for k in 1:K
	
		if (OS_NAME != :Windows) print("$k $V\r") end

 		DW = randn(2, N-1) .* sqrt(dt)
	 	yy = LinProc.eulerv(0.0, u, b, sigma, Dt, DW)


	 	
		v = yy[1:2, N]

		v0 = re()

		
 		# find lambda at the endpoint
		lambda = Lyap.lyap(B', -a(T,v0))
		lplambda = LinProc.lp(T, u, v0, B, beta, lambda)
		
 		#yy = LinProc.eulerv(0.0, u, v0, LinProc.Bcirc(T, v0, b, sigma, B, beta, lambda), sigma, Dt, DW)
	 	#ll = LinProc.llikelixcirc(0, T, yy, b, a, B, beta, lambda)
		DW = randn(2, N-1) .* sqrt(ds)
		u0 =  LinProc.UofX(0,u,  T, v0,  B, beta)
		yy = LinProc.eulerv(0.0, u0, LinProc.bU(T, v0, b, a, B, beta, lambda),  (s,x) -> sqrt(T)*sigma(LinProc.ddd(s,x, 0.0, T, v0,  B, beta)...), Ds, DW)
		
		ll = LinProc.llikeliU(S, yy, T, v0, b, a,  B, beta, lambda)
	
 

		l =  exp(lplambda + ll)/pe(v0)
		
 		# running mean and sum of squares
 		if (ll > llmax[2]) vmax = v0 end
		llmax = [min(llmax[1], ll), max(llmax[2], ll)]
		lplambdamax = [min(lplambdamax[1], lplambda), max(lplambdamax[2], lplambda)]
		 
		L += l
		L2 += l^2		
		V0 += E(v0) * l
		V02 += (E(v0) * l).^2
		V += 1.*E(v)
		V2 += E(v).^2
		#println("$L $V0 $V")    
		if (0 == k % 200)
			print("$k:")
		  
			p =  mc2(k, L, L2)
			println(" v ", mc2(k, V, V2)," v0 ", mc2(k,float64(V0), float64(V02)), " < p $p max's ll ", round(llmax, 1), " lp ", round(lplambdamax,1),", v0max ", round(vmax,2)," >" )
	 
		end	
	end
end

function dens(K, N, v0, t, T, B, A)

	println("Compute P(XT in dv), K=$K, T=$T, t=$t")
	println("E ~p(t, Xt)  =  ~p(0,u) * E exp(D(X°)(t)) ")
 
	llmax = [Inf,-Inf]
	llxmax = [Inf,-Inf]


	# avoid cancellation 
	Lx = BigFloat(0.0)
	Lx2 = BigFloat(0.0) 
	Lo = BigFloat(0.0)
	Lo2 = BigFloat(0.0) 

	
	Dt = diff(linspace(0., t, N))
	dt = Dt[1]
	Smax = LinProc.taui(t,T)
	S = linspace(0.,Smax, N)
	Ds = diff(S)
	ds = Ds[1]
	lambda = Lyap.lyap(B', -A)
	lplambda = LinProc.lp(T, u, v0, B, beta, lambda)
		
	for k in 1:K
	
		if (OS_NAME != :Windows) print("$k \r") end

		DW = randn(2, N-1) .* sqrt(dt)
	 	yy = LinProc.eulerv(0.0, u, b, sigma, Dt, DW)
		v = yy[1:2, N]
		llx = LinProc.lp(T-t, v, v0, B, beta, lambda)
		

		 
 		DW = randn(2, N-1) .* sqrt(ds)
		u0 =  LinProc.UofX(0,u,  T, v0,  B, beta)
		yy = LinProc.eulerv(0.0, u0, LinProc.bU(T, v0, b, a, B, beta, lambda),  (s,x) -> sqrt(T)*sigma(LinProc.ddd(s,x, 0.0, T, v0,  B, beta)...),  Ds, DW)
		
		ll = LinProc.llikeliU(S, yy, T, v0, b, a,  B, beta, lambda)
	
		lo =  exp(lplambda + ll)
		
		# running mean and sum of squares
		llmax = [min(llmax[1], ll), max(llmax[2], ll)]
		llxmax = [min(llxmax[1], llx), max(llxmax[2], llx)]

		Lo += lo
		Lo2 += lo^2		
		Lx += exp(llx)
		Lx2 += exp(2*llx)
		#println("$L $V0 $V")    
		if (0 == k % 1)
			print("$k:")
		  
			println(" ", mc2(k, float64(Lx), float64(Lx2)), " ~ ", mc3(k,float64(Lo), float64(Lo2)), " < max's llo ", round(llmax, 1), " llx ", round(llxmax,1)," >" )
	 
		end	
	end
end




N = 601 #design points
K = 1E6 #samples

#v = 1.15 .* LinProc.mu(T, u, B, beta)
v = 1.4*LinProc.mu(T, u, B, beta)
#v = [0.35373318326262687,-0.15995621287918171]
#test, whether X(T) in E
ep = 0.05
function E(x)
 max(abs(x - v)) <= ep
end

# a density, supported on E
function pe(x)
 E(x) / (2*ep)^2
end

# sampling from pe
function re()
  v + 2*ep*rand(2)-ep
end

#rare(K, N, E, re, pe)
println("D()=\ndens(K, N, v, 0.9999*T, T, B, 1.0*a(T,v))")

	D() = dens(K, N, v, 0.999*T, T, B, 1.*a(T,v))
#K = 1E4
#som = 0.
#dnorm(x) = exp(-0.5*x*x)/sqrt(2pi)
#dkreis(x) = (norm(x)< 1)*real(1./pi)

#for i in 1:K
#	x = randn(2,1)
#	som += pe(x)/dnorm(x[1])/dnorm(x[2])

#end
#print(som/K)

#println(exp(LinProc.lp(T, u, v, B, beta, lambda0)))


