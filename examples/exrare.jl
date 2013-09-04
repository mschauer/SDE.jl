using Lyap

require("linproc.jl")
require("misc.jl")

srand(3)

 
th = 1.3
si = 0.1
# si = 0.05
u = [-sqrt(th), -0.5] # start below focus
d = 2


b(s,x) = [x[2], -x[2] + (th - x[1]*x[1])*x[1]]
#b has root/focus at [-sqrt(th), 0]
#linearization of (th - x[1]*x[1])*x[1] in -sqrt(th): y  =  -2*th*x[1]-2*th^(3/2)
B = [0 1; -2*th -1]
beta = [0, -2*th^(3/2)]

function sigma(s,y)
	x = copy(y) - [-sqrt(th),0]
	m = norm(x)
	if (m <= eps()) return(sqrt(2.0) * 0.5 .* [[1. -1.],[1. 1.]]) end
	rho = 1 + 5*atan(m)
	si/m*[[x[2], -x[1]]  [rho*x[1], rho*x[2]]]
end
#sigma(s,x) = si*eye(2)

a= (s,x) -> sigma(s,x)*sigma(s,x)'



# if v=e() is distributed with log density le()
# and yy is a path with drift LinProc.Bcirc(T, v, b, sigma, B, beta, lambda) and diffusion sigma
# (lambda has to be computed for B(T, v) and a(T,v)!)
# then the expectation of llikelicirc(... yy, ...) is 1/pbar
# where pbar = exp(lp(..., B, beta, lambda) - le(v0)))
		
#7900: p [0.993,0.066] max's [-23.3,21.7][-48.0,5.3] v [-1.22742,-0.38162,0.00113,0.00266] v0 [-1.222,-0.39,0.0794,0.0202]

#T = 0.2, N = 1201, K = 1E6
#200000: v [0.02534,0.00069] v0 [0.02702,0.00034] < p [0.02702,0.00034] max's ll [-75.6,4.4] lp [-0.5,2.8] >
function rare(K, N, E, e, pe)

	println("Compute P(X in E), K=$K")
 
	llmax = [Inf,-Inf]
	lpbarmax = [Inf,-Inf]

	L = L2 = 0.0
	V = 0.0
	V2 = 0.0

	# avoid cancellation 
	V0 = BigFloat(0.0)
	V02 = BigFloat(0.0)
	
	Dt = diff(linspace(0., T, N))
	dt = Dt[1]
	
	for k in 1:K
	
		if (OS_NAME != :Windows) print("$k $V\r") end
		#sample endpoint of tilde X with lambda0
#		v0 = LinProc.sample_p(T, u, B, beta, lambda0) 

		DW = randn(2, N-1) .* sqrt(dt)
	 	yy = LinProc.eulerv(0.0, u, b, sigma, Dt, DW)
		v = yy[1:2, N]

		v0 = e()

		
 		# find lambda at the endpoint
		lambda = Lyap.lyap(B', -a(T,v0))
		lpbar = LinProc.lp(T, u, v0, B, beta, lambda)
		
 		DW = randn(2, N-1) .* sqrt(dt)
	 	yy = LinProc.eulerv(0.0, u, v0, LinProc.Bcirc(T, v0, b, sigma, B, beta, lambda), sigma, Dt, DW)
	 	ll = LinProc.llikelixcirc(0, T, yy, b, a, B, beta, lambda)
		l =  exp(lpbar + ll)/pe(v0)
		#println(pbar)
		#print(v0, yy[1:10], ll)
		# running mean and sum of squares
		llmax = [min(llmax[1], ll), max(llmax[2], ll)]
		lpbarmax = [min(lpbarmax[1], lpbar), max(lpbarmax[2], lpbar)]

		L += l
		L2 += l^2		
		V0 += E(v0) * l
		V02 += (E(v0) * l).^2
		V += 1.*E(v)
		V2 += E(v).^2
		#println("$L $V0 $V")    
		if (0 == k % 100)
			print("$k:")
		  
			p =  mc2(k, L, L2)
			println(" v ", mc2(k, V, V2)," v0 ", mc2(k,float64(V0), float64(V02)), " < p $p max's ll ", round(llmax, 1), " lp ", round(lpbarmax,1)," >" )
	 
		end	
	end
end


K = 40000
T = 0.2


N = 1201 #design points
K = 1E6 #samples

#v = 1.15 .* LinProc.mu(T, u, B, beta)
v = 1.15*LinProc.mu(T, u, B, beta)

#test, whether X(T) in E
function E(x)
 max(abs(x - v)) <= 0.05
end

# a density, supported on E
function pe(x)
 E(x) * 100.0
end

# sampling from pe
function e()
  v + 0.10*rand(2)-0.05
end

rare(K, N, E, e, pe)


#println(exp(LinProc.lp(T, u, v, B, beta, lambda0)))


