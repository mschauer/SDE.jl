module ExBridge
using Winston
using SDE
using LinProc
include(Pkg.dir("SDE","src", "leading.jl"))
include(Pkg.dir("SDE","src", "misc.jl"))

d = 2 
SV = false #save images?
PLOT = true
si = 0.1; include("excoeff2.jl")
K = 100
example = ["brown", "lin"][1]
letters =  map(string,char(int('A'):int('Z')))
#plot proposals in all segments

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
 
	p
end

#plot proposals in one segment
function plstep(p, xt, xd, yprop, m, nh, letters)
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

 	xprop = yprop[m]
	
	add(p, Curve(xprop[1,:],xprop[2,:], "color","blue", "linewidth", 0.1))
 
	add(p, Points(xd[1,m:(m+1)],xd[2,m:(m+1)],"symboltype", "filled circle", "color", "red", "size", 1.2))
	add(p, Points(xd[1,m:(m+1)],xd[2,m:(m+1)],"symboltype", "empty circle", "color", "white", "size", 0.6))

	add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "white", "size", 2.5))
	add(p, Winston.DataLabel(1.1*xd[1,m+1],1.1*xd[2,m+1]+0.1,  letters[m+1], "color", "white", "size", 2.5))
	add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "white", "size", 1.5))
	add(p, Winston.DataLabel(1.1*xd[1,m+1],1.1*xd[2,m+1]+0.1,  letters[m+1], "color", "white", "size", 1.5))
	add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "black", "size", 2))
	add(p, Winston.DataLabel(1.1*xd[1,m+1] ,1.1*xd[2,m+1]+0.1,  letters[m+1], "color", "black", "size", 2))
	

	#add(p, Points(xprop[1,nh],xprop[2,nh],"symboltype", "filled circle", "color", "red", "size",  0.7*min(5.,(ll[m]/llmax[m])))) #midpoints
	#add(p, Curve(xt[1,:],xt[2,:], "color","grey")) #true function
	 
	p
end

#plot proposals weighted
function plstep(p, xt, xd, yprop, ll, llmax, m, nh, letters)
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

 	xprop = yprop[m]
	
	#add(p, Curve(xt[1,:],xt[2,:], "color","grey"))
	add(p, Curve(xprop[1,:],xprop[2,:], "color","blue", "linewidth", 0.2))
	add(p, Curve(xprop[1,:],xprop[2,:], "color","red", "linewidth", 0.7*min(5.,(ll[m]/llmax[m]))))
#	add(p, Curve(xprop[1,:],xprop[2,:], "color","blue", "linewidth", 0.7*min(5.,(1.-ll[m]/llmax[m]))))
	 
	add(p, Points(xd[1,m:(m+1)],xd[2,m:(m+1)],"symboltype", "filled circle", "color", "red", "size", 1.2))
	add(p, Points(xd[1,m:(m+1)],xd[2,m:(m+1)],"symboltype", "empty circle", "color", "white", "size", 0.6))

	add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "white", "size", 2.5))
	add(p, Winston.DataLabel(1.1*xd[1,m+1],1.1*xd[2,m+1]+0.1,  letters[m+1], "color", "white", "size", 2.5))
	add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "white", "size", 1.5))
	add(p, Winston.DataLabel(1.1*xd[1,m+1],1.1*xd[2,m+1]+0.1,  letters[m+1], "color", "white", "size", 1.5))
	add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "black", "size", 2))
	add(p, Winston.DataLabel(1.1*xd[1,m+1] ,1.1*xd[2,m+1]+0.1,  letters[m+1], "color", "black", "size", 2))
	
	#add(p, Points(xprop[1,nh],xprop[2,nh],"symboltype", "filled circle", "color", "red", "size",  0.7*min(5.,(ll[m]/llmax[m])))) #midpoints
	 
	p
end
function plobs(xt, xd, letters)
	p = FramedPlot()
	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)
	
	add(p, Curve(xt[1,:],xt[2,:], "color","black", "linewidth", 0.5))
	add(p, Points(xd[1,:],xd[2,:],"symboltype", "filled circle", "color", "red", "size", 1.2))
	add(p, Points(xd[1,:],xd[2,:],"symboltype", "empty circle", "color", "white", "size", 0.6))

	for m in 1:length(letters)
		add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "white", "size", 2.5))
		add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "white", "size", 1.5))
		add(p, Winston.DataLabel(1.1*xd[1,m],1.1*xd[2,m]+0.1, letters[m], "color", "black", "size", 2))
	end	 
	
	p
end
 



 
############ 
print("Parameters: ")
############ 

M = 3 #number of bridges
n = 1000 #samples each bridge
N = 12000 + 1   #full observations
TT = 4*M		#time span
mar = 0.3
#example = "lin"

println("M$M n$n N$N T$TT $example")

############

############
print("Generate X")
############

srand(6)
Dt = diff(linspace(0., TT, N))
dt = Dt[1]
DW = randn(2, N-1) .* sqrt(dt)
x = LinProc.eulerv(0.0, u, b, sigma, Dt, DW)

	

#compute range of observations
R1 = range(x[1,:]) 

R1 = (R1[1] -mar*(R1[2]-R1[1]),R1[2] + mar*(R1[2]-R1[1]))
R2 = range(x[2,:]) 
R2 = (R2[1]  -mar*(R2[2]-R2[1]), R2[2] + mar*(R2[2]-R2[1]))


#R1 = (-10., 10.)
#R2 = (-10., 10.)  
println(".")

xd = x[:, 1:(N-1)/M:end]
xtrue = x[:, 1:(N-1)/M/n:end]
if PLOT
	p = plobs(x, xd, ["A", "B", "C", "D"])
	Winston.display(p)
	if(SV) Winston.file(p, "img/obs.png") end
end

 

############
println("Generate bridges (using proposals $example).")
############
srand(3)
 
yprop = cell(M)
y = cell(M)
U = zeros(2,n+1)
llold = zeros(M)
ll = zeros(M)
llmax = zeros(M)
yprop = cell(M)
y = cell(M)
yy = zeros(2,n+1)

 
 S = Smax =  Y = Tmax = Tmin = []
 
#accepted bridges
bb = zeros(M)


p = cell(M)
pprop = cell(M)
if PLOT 
	for m = 1:M 
		p[m] = FramedPlot()
		pprop[m] = FramedPlot()
	end
end

for k = 1:K
 	dt = TT/M/n
	if (k == K/2 && SV)
		srand(3)
	end
	
	for m = 1:M
		u = leading(xd, m)
		v = leading(xd, m+1)
		
		Tmax = m/M*TT 
		Tmin = (m-1)/M*TT
		if example == "brown"; 
			B = 0*bB;  
			beta = 0*bbeta
		end
		if example == "lin"; B = exp(-0.2*Tmax)*bB; end
		
		T = Tmax - Tmin
		Ts =  linspace(0, T-dt, n)
		lambda = inv(a(Tmax,v))

		try
			lambda = Lyap.lyap(B', -a(Tmax,v))
		end
	 	Smax = LinProc.taui(T-2*dt,T)
		S = linspace(0.0,Smax, n)

	
		#S = -log(1.-Ts/T) 
		Ds = diff(S)
		nh = int(LinProc.taui(0.5T,T)/Ds[1])

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

		if PLOT
			if (k > K/2  )
				p[m] = plstep(p[m], xtrue, xd, yprop, ll, llmax,  m, nh, letters)
			else	
				llmax = max(llmax, ll) 
				pprop[m] = plstep(pprop[m], xtrue, xd, yprop,  m, nh, letters)
			end
		else 
			llmax = max(llmax, ll) 
		end

	end
	println("\t",round(llmax,3))
	if PLOT
		if (k > K/2)
	 		Winston.display(p[1])
		else
		 	Winston.display(pprop[1])
		end		
	end


	println("k $k\t acc ",  round(100*mean(bb)/k,2), round(100*bb/k,2), )

end
if(SV)
for m = 1:M; Winston.file(p[m], "img/$example$m.png") end 
for m = 1:M; Winston.file(pprop[m], "img/$(example)prop$m.png") end 
end 


println("bridge acc %",  round(100*mean(bb)/K,2), round(100*bb/K,2), )
end
