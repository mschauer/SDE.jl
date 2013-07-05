using Lyap
#using Diffusion
require("leading.jl")
require("linproc.jl")
require("misc.jl")
srand(3)

SV = false #save images?


function ito(x, dy)
	n = length(dy) + 1
	y = 0.0
	for i in 2:n
		y = y + x[i-1]*dy[i-1] 
	end
	y
end
  

function euler(t0, u, b, sigma, dt, dw::Matrix)
	S = size(dw)
	endd = length(S)	
	N = S[end] + 1
	
	shape = size(sigma(0,u)*leading(dw,1))
 
	X = zeros(shape..., N)

	#delta = copy(u)
	y = copy(u)
	t = t0

	
	for i in 1:N-1
		subleading(X,i)[:] = y
		t += dt[i]
		y[:] = y .+  b(t,y)*(dt[i]) .+ sigma(t,y)*leading(dw, i)
	
	end
	subleading(X,N)[:] = y
	X
end
function eulerv(t0, u, v, b, sigma, dt, dw::Matrix)
	X = euler(t0, u,  b, sigma, dt, dw)
	X[:,end] = v
	X
end

function likelixcirc(t, T, v, Xcirc, b, a,  B, beta, lambda)
#ignores last value
	
	function L(s,x)
		R = LinProc.H(T-s, B, lambda)*(x - LinProc.V(T-s, v, B, beta))
	  	return (b(s,x) - B*x - beta)' * R + 0.5 *trace((a(s,x) - a(T,v)) *( LinProc.H(T-s, B, lambda) + R*R'))
	end
	
	sum = 0
	N = size(Xcirc,2)
	s= t
	x = v 
	for i in 0:N-1-1 #skip last value, summing over n-1 elements
	  s = t + (T-t)*(i)/(N-1) 
	  x = leading(Xcirc, i+1)
	  sum += scalar(L(s, x)) * (T-t)/N
	end
	sum += scalar( 2*sqrt(T-s)*(b(T,x) - B*x - beta)'* a(T,v)*(v-x)) #interpolate drift part of last interval like square root
	
	exp(sum)
end

function pl(x)
	p = FramedPlot()
	setattr(p, "xrange", (-1.5,-1))
	setattr(p, "yrange", (-0.5,0.5))

	add(p, Curve(x[1,:],x[2,:]))
	Winston.display(p)

	p	
end

function plstep(xt, xd, y, yprop)
	p = FramedPlot()
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

	x = apply(hcat, y)
	xprop = apply(hcat, yprop)
	
	add(p, Curve(xt[1,:],xt[2,:], "color","grey"))
	add(p, Curve(xprop[1,:],xprop[2,:], "color","light blue"))
	add(p, Curve(x[1,:],x[2,:], "color","black"))

	add(p, Points(xd[1,:],xd[2,:],"type", "dot", "color","red"))
	
	Winston.display(p)
	 
	p
end

function plobs(xt, xd)
	p = FramedPlot()
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)
	
	add(p, Curve(xt[1,:],xt[2,:], "color","black"))
	add(p, Points(xd[1,:],xd[2,:],"type", "filled circle", "color","red"))
	
	Winston.display(p)
	 
	p
end

#MC estimate with normal 95% confidence assuming independent samples

 

#mc(Z) == mc2(length(Z),sum(Z), sum(Z.^2))  

#th = 1.7
#b has root/focus at [-sqrt(th), 0]
#linearization of (th - x[1]*x[1])*x[1] in -sqrt(th): y  =  -2*th*x[1]-2*th^(3/2)
th = 1.3
si = 0.1
# si = 0.05
u = [-sqrt(th), -0.5] # start below focus
d = 2
zd = zeros(d)
od = ones(d)
Id = eye(d)

b(s,x) = [x[2], -x[2] + (th - x[1]*x[1])*x[1]]
#Lb(s, x) = [x[2], -x[2] - 2*th*x[1] - 2*th^(3/2)] 

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

K = 40000
println("Compute p(x,y), K=$K")
T = 0.8
v = [-1.3, 0]
tb(s, x) = B*x + beta
tsigma(s,x) = sigma(T, v)
ta(s,x) = tsigma(s, x)*tsigma(s, x)'
lambda = Lyap.lyap(B', -a(T,v))
 
N = 1001 #samples
Dt = diff(linspace(0., T, N))
dt = Dt[1]

function kern(z, d, h)
         (exp (-1/2*d*log(2pi*h)  -0.5*z.*z/h))
         
end
lambda0 = Lyap.lyap(B', -a(T,u))
 		

LL = BigFloat(0.)
LL2 = BigFloat(0.)

# if the likehood function is correct, then this MC estimate converges to 1.
# if v_0 is sample_p(T, u, B, beta, lambda0) distributed
# and yy is a path with drift LinProc.Bcirc(T, v0, b, sigma, B, beta, lambda) and diffusion sigma
# then the expectation of likelicirc(... yy, ...) is 1/pbar
# where pbar = exp(lp(..., B, beta, lambda) - lp(... B, beta, lambda0)))
		

function testlikeli1(K, N, T)
 
	LL = LL2 = 0.0
	Dt = diff(linspace(0., T, N))
	dt = Dt[1]

	for k in 1:K
	
		if (OS_NAME != :Windows) print("$k\r") end
		#sample endpoint of tilde X with lambda0
		v0 = LinProc.sample_p(T, u, B, beta, lambda0) 
		
 		# find lambda at the endpoint
		lambda = Lyap.lyap(B', -a(T,v0))
		# correction for change in lambda
		pbar = scalar(exp(LinProc.lp(T, u, v0, B, beta, lambda)-LinProc.lp(T, u, v0, B, beta, lambda0)))
		
 		DW = randn(2, N-1) .* sqrt(dt)
	 	yy = eulerv(0.0, u, v0, LinProc.Bcirc(T, v0, b, sigma, B, beta, lambda), sigma, Dt, DW)
		ll =  pbar*scalar(likelixcirc(0, T, v0, yy, b, a, B, beta, lambda))
		#println(pbar)
		#print(v0, yy[1:10], ll)
		# running mean and sum of squares
		LL += ll
		LL2 += ll^2		
		#println("$k $v0 $ll $LL $LL2")    
		if (0 == k % 100)
			print("$k:")
		  
			p =  mc2(k, LL, LL2)
			println(": p \t $p")
	 
		end	
	end
end


#if LinProc.lp() is correct, then this MC estimate converges to 1.

function testlikeli2(K, N, T)


	LL = LL2 = 0.0
	Dt = diff(linspace(0., T, N))
	dt = Dt[1]
	mu =  B*u+beta
	gamma = inv(a(T,v))


	for k in 1:K
		if (OS_NAME != :Windows) print("$k\r") end
 	 	DW = randn(2, N-1) .* sqrt(dt)
	 	yy = euler(0.0, u, tb, tsigma, Dt, DW)
		y = yy[1:2, N]
		
		ll = scalar(exp(LinProc.lp0(T, u, y, mu, gamma)-LinProc.lp(T, u, y, B, beta, lambda)))

		# running mean and sum of squares
		LL += ll
		LL2 += ll^2	

  		if (0 == k % 100)
			print("$k:")
		  
			p =  mc2(k, LL, LL2)
			println(": p \t $p")
	 
		end	
	end
end


N = 201 #design points
K = 1E6 #samples
testlikeli1(K, N, 0.2)

