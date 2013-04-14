
#using Diffusion
include("leading.jl")



function ito(x, dw)
	n = length(dw) + 1
	y = 0.0
	for i in 2:n
		y = y + x[i-1]*dw[i-1] 
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


th = 1.98
#th = 1.3
u = [-1., 0.]

b(s,x) = [x[2], -x[2] + (th - x[1]*x[1])*x[1]]
si = 0.02
#si = 1.
sigma(s,x) = diagm([si, si])


N = int(1e4)
T = 10 #int(1e2)
dt = diff(linspace(0., T, N))
dw = randn(2, N-1) .* sqrt(dt[1])
#euler(t0, u, b, sigma, dt::Float64, dw::Vector)
x = euler(0.0, u, b, sigma, dt, dw)

#Prior variance
s = 1


function plotx(x)
	p = FramedPlot()
	add(p, Curve(x[1,:],x[2,:]))
	Winston.display(p)
	
end


L(theta, x) = ito(-x[2, :] + theta*x[1, :] - x[1, :].^3, diff(x[2, :],2))- 0.5*ito( (-x[2, :] + theta*x[1, :] - x[1, :].^3).^2 , dt)

m(x) = ito(x[1, :], diff(x[2, :],2))/(si^2) - ito(x[1, :].*(- x[2, :] - x[1, :].^3), dt)/(si^2)
mu = m(x)


W =  ito(x[1, :].^2, dt)/(si^2)  +  1./s^2
#mu = ito(x[1, :], diff(x[2, :],2))/(si^2)  + ito(x[1, :] .* x[2, :] .+ x[2, :].^4, dt )/(si^2) wrong!

println("th= ", mu/W, "+-", 1/W)

