using Lyap
#using Diffusion
require("leading.jl")
require("linproc.jl")
srand(3)
range(x) = (min(x), max(x))

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

function scalar(x)
	assert(length(x) == 1)
	x[1]
end

function likelixcirc(t, T, v, Xcirc, b, a,  B, beta, lambda)
	
	function L(s,x)
		R = LinProc.H(T-s, B, lambda)*(x - LinProc.V(T-s, v, B, beta))
	  	return (b(s,x) - B*x - beta)' * R + 0.5 *trace((a(s,x) - a(T,v)) *( LinProc.H(T-s, B, lambda) + R*R'))
	end
	
	sum = 0
	N = size(Xcirc,2)
	for i in 1:N
	  s = t + (T-t)*(i-1)/N
	  x = leading(Xcirc, i)
	  sum += scalar(L(s, x)) * (T-t)/N
	end
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

function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, 1.96*ste], int(2-log(10,1.96*ste)))
end


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


K = 4E4
println("Compute p(x,y) ($K iterations)")
T = 0.8
v = [-1.3, 0]
tb(s, x) = B*x + beta
tsigma(s,x) = sigma(T, v)
ta(s,x) = tsigma(s, x)*tsigma(s, x)'
lambda = Lyap.lyap(B', -a(T,v))
 
N = 1001 #samples
Dt = diff(linspace(0., T, N))
dt = Dt[1]
Z = [1]
Z0 = [1]
LL = [1]
for k in 1:K
	print("$k\r")
	DW = randn(2, N-1) .* sqrt(dt)
	x = norm(v-leading(euler(0.0, u,  b,  sigma, Dt, DW), N))
	#z = scalar(exp(LinProc.lp0(1E-5,  z, v, zd, Id)))
	Z = vcat(Z,x) 
	x = norm(v-leading(euler(0.0, u,  tb,  tsigma, Dt, DW), N))
 	Z0 = vcat(Z0,x)
	if(k < K/4) 
	
	 	yy = euler(0.0, u, LinProc.Bcirc(T, v, b, sigma, B, beta, lambda), sigma, Dt, DW)
		ll =  scalar(likelixcirc(0, T, v, yy, b, a, B, beta, lambda))
		LL = vcat(LL,ll)
	end
		            
	#println("x(T) = $(leading(x, N)), kernel(x(T)) = $z")
end
function kern(z, d, h)
         (exp (-1/2*d*log(2pi*h)  -0.5*z.*z/h))
         
end
#@vectorize_1arg Number kern

#function kernel_est(Z, h)
# M = map(z -> kern(2,h, z),Z)
# mc(M)
#end

h = 1/K^0.9
pnaiv =  mc(kern (Z ,d, h))
p0naiv = mc(kern (Z0,d, h))
p0 = scalar(exp(LinProc.lp(T, u, v, B, beta, lambda)))
p =  mc(p0*LL)
 
#mean((Z .< 0.05)/0.05)/ mean((Z0 .< 0.05)/0.05)

println("h = $h\npnaiv\t$pnaiv \np0naiv\t$p0naiv \np\t$p \np0\t$p0")
println("LL \t", mc(LL), "\t\tpnaiv[1]/p0naiv[1]\t", pnaiv[1]/p0naiv[1])


for h in (logspace(log(10,1/sqrt(K)),log(10, 1/K^1.2), 10))
	pnaiv =  mc(kern(Z ,d, h))	;p0naiv = mc(kern(Z0 ,d, h))
	println("h = $h\npnaiv\t$pnaiv \np0naiv\t$p0naiv")
end
stop()


print("Generate x")

N = 12000 + 1   #full observations
T = 5		#time span
T = 4
Dt = diff(linspace(0., T, N))
dt = Dt[1]
DW = randn(2, N-1) .* sqrt(dt)
x = euler(0.0, u, b, sigma, Dt, DW)


#compute range of observations
R1 = range(x[1,:]) 
 
R1 = (R1[1] -0.1*(R1[2]-R1[1]),R1[2] + 0.1*(R1[2]-R1[1]))

R2 = range(x[2,:])
R2 = (R2[1]  -0.1*(R2[2]-R2[1]), R2[2] + 0.1*(R2[2]-R2[1]))


#R2 = (-1.5,-1)
#R2 = (-0.5,0.5)
println(".")


#Prior variance

s = 20.0



L(theta, x, Dt) = ito(-x[2, :] + theta*x[1, :] - x[1, :].^3, diff(x[2, :],2))- 0.5*ito( (-x[2, :] + theta*x[1, :] - x[1, :].^3).^2 , Dt)

#mu(x) = ito(x[1, :], diff(x[2, :],2))/(si^2) - ito(x[1, :].*(- x[2, :] - x[1, :].^3), dt)/(si^2)
mu0(x, Dt) = ito(x[1, :], diff(x[2, :],2))/(si^2) + ito(x[1, :].* x[2, :] + abs(x[1, :]).^4, Dt)/(si^2)
W0(x, Dt) =  ito(x[1, :].^2, Dt)/(si^2) 



m = mu0(x, Dt)
w = W0(x, Dt)  +  1./s^2
println("th= ", round(m/w,5), "+-", round(1.96*sqrt(1/w), 5))
 
###
M = 5 #number of bridges
n = 600 #samples each bridge
xd = x[:, 1:(N-1)/M:end]
xtrue = x[:, 1:(N-1)/M/n:end]

p = plobs(x, xd)
if(SV) Winston.file(p, "img/s$(si)obs.png") end

a=(s,x) -> sigma(s,x)*sigma(s,x)'

v = [-sqrt(th), 0.0]
B = [0 1; -2*th -1]
beta = [0, -2*th^(3/2)]
lambda = Lyap.lyap(B', -a(0,v))
z = eulerv(0.0, u, v,  LinProc.Bcirc(T, v, b, sigma, B, beta, lambda), sigma, Dt, DW)
ll = likelixcirc(0.0, T, v, z, b, a, B, beta, lambda)
#println(ll)
#error()

#th = 3/4*th
yprop = cell(M)
y = cell(M)
yy = zeros(2,n+1)
K = 100
thetas = zeros(K)
alpha = 1.0 


#accepted bridges
bb = 0
for k = 1:K
	b(s,x) = [x[2], -x[2] + (th - x[1]*x[1])*x[1]]
	th2 = th
	B = alpha*[0 1; -2*th2 -1]
	beta = alpha*[0, -2*th2^(3/2)]
	dt = T/M/n
	Dt = dt*ones(n-1)
	
	for m = 1:M
		u = leading(xd, m)
		v = leading(xd, m+1)
	#	println("u -> v $u $v")
		DW = randn(2,n-1)*sqrt(dt)
	
		lambda = Lyap.lyap(B', -a(0,v))

		
		yy = euler(((m-1)/M)*T, u, LinProc.Bcirc((m/M)*T, v, b, sigma, B, beta, lambda), sigma, Dt, DW)
		yprop[m] = yy
		if(k == 1) y[m] = yy end 
		ll =  likelixcirc(((m-1)/M)*T, (m/M)*T, v, yy, b, a, B, beta, lambda)
		llold = likelixcirc(((m-1)/M)*T, (m/M)*T, v, y[m], b, a, B, beta, lambda)
		#
		#println(ll)
		#readline(STDIN)
		if (llold > 0)
			acc = min(1.0, ll/llold)
		else
			acc = 1.0
		end
		println("\t acc ", round(acc, 3), " ", round(llold,3), " ", round(ll,3), " ", round(acc,3))
	
		if (rand() <= acc)
			  bb += 1
			  y[m] = yy
			 
		end	
			

	end
	p = plstep(xtrue, xd, y, yprop )
	if(SV && k > K/2) Winston.file(p, "img/s$(si)a$(alpha)k$k.png") end
	mu = 0.
	W = 1./s^2
		
	for m = 1: M
		mu += mu0(y[m], Dt)   
		W  += W0(y[m], Dt)     
	end
	#
	#th = mu/W + norm(1)*sqrt(1/W)	
	th=th
	thetas[k] = th

	println("k $k\tth $th acc ", round(100*bb/k/M,2 ))
end

thetas2 = thetas[(K/2):end]
mi = round(mean(thetas2), 3)
ci = round(1.96*std(thetas2)/sqrt(length(thetas2)),3)
println("theta = [", mi-ci, ", ", mi + ci, "] (95%-ci)")

println("bridge acc %", round(100*bb/K/M,2))

