#using Winston
 
#using Diffusion
require("leading.jl")
require("linproc.jl")
require("misc.jl")
require("lyap.jl")
srand(3)
 
SV = false #save images?



si = 0.08
d = 2
B = [ -0.2   -1;   0.5  -0.4]
bB = copy(B)
#B = [ -2.21225   -0.49590;   0.631753  -1.34795]
#A = [ 0.170252  0.178533; 0.178533  2.43255]

beta = [0.5,-0.5]
#beta = [0.,0.]

function sigma0(s,y)
	x = copy(y) - beta
	m = norm(x)
	if (m > eps())  
		rho = 1 + 5*2atan(m)/pi
		return( si*(eye(2)+ 2atan(m)/pi/m*[[x[2], -x[1]]  [rho*x[1], rho*x[2]]]))
	end
	return si*eye(2)
end
#sigma(s,x) = [[1., -1.]  [5., 5.]]
sigma(t,x) = sigma0(0., [sin(t), cos(t)])
#sigma(t,x) = sigma0(t,x)

a(s,x) = sigma(s,x)*sigma(s,x)'

 

b(t,x) = exp(-0.2*t)*bB*x + beta
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
#b(t,x) =  atan(norm(x))*B*x + beta
 
u = [1.,0.]
 
 


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
function plstep(p, xt, xd, yprop, ll, llmax, m)
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

 	xprop = yprop[m]
	
	#add(p, Curve(xt[1,:],xt[2,:], "color","grey"))
	add(p, Curve(xprop[1,:],xprop[2,:], "color","black", "linewidth", 0.7*min(5.,(ll[m]/llmax[m]))))
	 
	add(p, Points(xd[1,m:(m+1)],xd[2,m:(m+1)],"symboltype", "filled circle", "color", "red", "size", 1.2))
	add(p, Points(xd[1,m:(m+1)],xd[2,m:(m+1)],"symboltype", "empty circle", "color", "white", "size", 0.6))
	 
	p
end
function plobs(xt, xd)
	p = FramedPlot()
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)
	
	add(p, Curve(xt[1,:],xt[2,:], "color","black", "linewidth", 0.5))
	add(p, Points(xd[1,:],xd[2,:],"symboltype", "filled circle", "color", "red", "size", 1.2))
	add(p, Points(xd[1,:],xd[2,:],"symboltype", "empty circle", "color", "white", "size", 0.6))
	
	Winston.display(p)
	 
	p
end
 


print("Generate x")
 
###
M = 4 #number of bridges
n = 1000 #samples each bridge
N = 12000 + 1   #full observations
TT = 0.5*M		#time span
 
Dt = diff(linspace(0., TT, N))
dt = Dt[1]
DW = randn(2, N-1) .* sqrt(dt)
x = LinProc.eulerv(0.0, u, b, sigma, Dt, DW)

	

#compute range of observations
R1 = range(x[1,:]) 
mar = 0.5
R1 = (R1[1] -mar*(R1[2]-R1[1]),R1[2] + mar*(R1[2]-R1[1]))
R2 = range(x[2,:]) 
R2 = (R2[1]  -mar*(R2[2]-R2[1]), R2[2] + mar*(R2[2]-R2[1]))


#R1 = (-10., 10.)
#R2 = (-10., 10.)  
println(".")


#Prior variance
 

 


xd = x[:, 1:(N-1)/M:end]
xtrue = x[:, 1:(N-1)/M/n:end]

p = plobs(x, xd)
if(SV) Winston.file(p, "img/obs.png") end

#      [1869.16,515.23,170074.054,301.773]
#k 100    acc 25.0[44.0,26.0,23.0,7.0]
#bridge acc %25.0[44.0,26.0,23.0,7.0]
#        [1.976,1.849,1.49,6.175]
#k 100    acc 75.5[88.0,78.0,77.0,59.0]
#bridge acc %75.5[88.0,78.0,77.0,59.0]

#println(ll)
error()

#th = 3/4*th
yprop = cell(M)
y = cell(M)
U = zeros(2,n+1)
llold = zeros(M)
ll = zeros(M)
llmax = zeros(M)
yprop = cell(M)
y = cell(M)
yy = zeros(2,n+1)
K = 100
 
 S = Smax =  Y = Tmax = Tmin = []
#accepted bridges
bb = zeros(M)
example = "lin"
#example = "brown"

p = cell(M)	
for m = 1:M 
	p[m] = FramedPlot()
end

srand(3)
for k = 1:K
 	dt = TT/M/n
	if (k == K/2)
		for m = 1:M; p[m] = FramedPlot(); end
		srand(3)
	end
	
	for m = 1:M
		u = leading(xd, m)
		v = leading(xd, m+1)
		
		Tmax = m/M*TT 
		Tmin = (m-1)/M*TT
		if example == "brown"; B = bB*0.01;  end
		if example == "lin"; B = exp(-0.2*Tmin)*bB; end
		
		T = Tmax - Tmin
		Ts =  linspace(0, T-dt, n)
		lambda = zero(a(Tmax,v))

		try
			lambda = Lyap.lyap(B', -a(Tmax,v))
		end
		#print(lambda)
		Smax = LinProc.taui(T-2*dt,T)
		S = linspace(0.0,Smax, n)

		#S = -log(1.-Ts/T) 
		Ds = diff(S)

	 	DW = randn(2,n-1).*[sqrt(Ds)  sqrt(Ds)]'
		
		u0 = LinProc.UofX(0.0,u,  T, v,  B, beta)
				
	 	U = LinProc.eulerv(0.0, u0, LinProc.bU(Tmin, T, v, b, a, B, beta, lambda), (s,x) -> sigma(LinProc.ddd(s,x, Tmin, T, v,  B, beta)...), Ds, sqrt(T)*DW)
	
		ll[m] = exp(LinProc.llikeliU(S, U, Tmin, T, v, b, a,  B, beta, lambda))
		time, Y = LinProc.XofU2(S, U,  Tmin, T, v,  B, beta) 
	 	yprop[m] = Y
	 	if(k == 1) 
			y[m] = Y
			llold[m] = ll[m]
		end 
		
		if (llold[m] > 0)
			acc = min(1.0, ll[m]/llold[m])
		else
			acc = 1.0
		end
		println("\t acc ", round(acc, 3), " ", round(llold[m],7), " ", round(ll[m],7), " ", round(acc,3))
	
		if (rand() <= acc)
			  bb[m] += 1
			  y[m] = Y
			  llold[m] = ll[m]
			 
		end	
		
		if (k < K/2) llmax = max(llmax, ll) end
	
		p[m] = plstep(p[m], xtrue, xd, yprop, ll, llmax,  m)
	

	end
	println("\t",round(llmax,3))

 	Winston.display(p[1])

		



	println("k $k\t acc ",  round(100*mean(bb)/k,2), round(100*bb/k,2), )

end
if(SV)
for m = 1:M; Winston.file(p[m], "img/$example$m.png") end 
end 


println("bridge acc %",  round(100*mean(bb)/K,2), round(100*bb/K,2), )
#Winston.display(p[1])
