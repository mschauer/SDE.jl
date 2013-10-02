using Winston
using SDE
using Diffusion
using LinProc
include(Pkg.dir("SDE","src", "misc.jl"))

pic = 1
function dp(p)
	Winston.display(p)
	Winston.file(p, "examples/img/talk$pic.png")
	global pic
	pic += 1
end

function lyap(B::Real, ma::Real)
    ma/(2.*B)
end


function rgb(r,g, b)
	r = cut(r,0.,1.)
	g = cut(g,0.,1.)
	b = cut(b,0.,1.)
 	prod(("#",
	repr(uint8(round(r*255.)))[3:],
	repr(uint8(round(g*255.)))[3:],
	repr(uint8(round(b*255.)))[3:]))
end


function addpath1(p, t, y)
	
	add(p, Curve(t,y, "color","black", "linewidth", 20./M))
	
	p
end

function addpath1(p, t, y, l)
	l = cut(l, 0.0, 1.0)
	add(p, Curve(t,y, "color",rgb(l, 0., 1.-l, ), "linewidth", 20./M*(l+1.)/2.))
	#print(rgb(l, 1.-l, 1. ), " ")
	p
end
#M = 2000
M = 100

T = 0.25
N = 1000

R1 = (0., T)
R2 = (-2., 3.)
ts = linspace(0., T,N)
dt = diff(ts)
B = -5.
sigma = 2.
a = sigma^2
beta = 0.
lambda = lyap(B, -a)
p = FramedPlot()
setattr(p, "xrange", R1)
setattr(p, "yrange", R2)

for m in	 1:M

	y = euler(0., 1., (t,x)-> B*x + beta, (t,x) -> sigma, dt)
	p = addpath1(p, ts, y)
end
dp(p)

srand(1)
p = FramedPlot()
setattr(p, "xrange", R1)
setattr(p, "yrange", R2)
for m in 1:M

	y = sigma*brown1(1./sigma, T, N)
	p = addpath1(p, ts, y)
end
dp(p)


p = FramedPlot()
setattr(p, "xrange", R1)
setattr(p, "yrange", R2)

srand(1)
llmax = -Inf
for m in 1:M

	y = 1. .+ sigma*brown1(0., T, N)
	ll = sum( 1/a*(B.*y[1:end-1] + beta).*diff(y)) + 0.5*T*a 
	llmax = max(llmax, ll)
	#print(round(exp(ll)), " ")
end
println(exp(llmax))
srand(1)
for m in 1:M

	y =  1. .+ sigma*brown1(0., T, N)
	ll = sum( 1/a*(B.*y[1:end-1] + beta).*diff(y)) + 0.5*T*a

	p = addpath1(p, ts, y, exp(ll-llmax))
end
dp(p)
 	
srand(1)
p = FramedPlot()
setattr(p, "xrange", R1)
setattr(p, "yrange", R2)
su0 = 0.
su0v = 0.
for m in 1:M

	y = sigma*bb(1./sigma, 0.5/sigma, T, N)
	su0 +=  y[N/2]
	su0v +=  y[N/2]^2
	p = addpath1(p, ts, y)
end
(z0, se0) = mc2(M, su0, su0v)
add( p, SymmetricErrorBarsY(T/2, [z0], [se0]) )
dp(p)


p = FramedPlot()
setattr(p, "xrange", R1)
setattr(p, "yrange", R2)
srand(1)
llmax = -Inf
lmean = 0.
for m in 1:M

	y = sigma*bb(1./sigma, 0.5/sigma, T, N)
	ll = sum( 1/a*(B.*y[1:end-1] + beta).*diff(y) - 0.5/a*(B.*y[1:end-1] + beta).^2 .* dt )
	#print(round(exp(ll)), " ")
	llmax = max(llmax, ll)
	lmean += exp(ll)
	
end
lmean /= M
println(exp(llmax))
srand(1)

Z = zeros(2)
Z2 = zeros(2,2)

for m in 1:M

	y = sigma*bb(1./sigma, 0.5/sigma, T, N)
	ll = sum( 1/a*(B.*y[1:end-1] + beta).*diff(y) - 0.5/a*(B.*y[1:end-1] + beta).^2 .* dt )
	
	Z += [y[N/2] * exp(ll), exp(ll)]
	Z2 += [y[N/2] * exp(ll), exp(ll)]*[y[N/2] * exp(ll) exp(ll)]
	p = addpath1(p, ts, y, exp(ll-llmax))
end
(z1, se1) = selfadjmc(M, Z, Z2)
add( p, SymmetricErrorBarsY(T/2, [z1], [se1]) )
dp(p)


srand(1)
p = FramedPlot()
setattr(p, "xrange", R1)
setattr(p, "yrange", R2)

su2 = 0.
su2v = 0.

for m in 1:M

	
	
	y = euler(0., 1., (s,x) -> scalar(LinProc.bstar(T, 0.5, B, beta, a, lambda)(s,x)), (t,x) -> sigma, dt)
	su2 += y[N/2]
	su2v += y[N/2]^2
	p = addpath1(p, ts, y)
	
	
end
(z2, se2) = mc2(M, su2, su2v)
add( p, SymmetricErrorBarsY(T/2, [z2], [se2]) )
dp(p)



println("Expect. at T/2: unweighted ($z0 $se0) weighted  ($z1 $se1) bridge  ($z2 $se2) ")
println("Average weight: ", lmean, " theoretical ", exp(lp(T, 1.0, 0.5, B, beta, lambda*eye(1)))/exp(LinProc.lp0(T, 1.0, 0.5, 0, 1/a*eye(1))))
