module Diffusion
# lines starting with "#%" contain the documentation in ReST
# the leading "#%" and up to 2 spaces are removed. 
# careful, depending on your tab length, 8 spaces and one tab may look identically aligned, but may end up not so
# in the corresponding style file



#%  .. currentmodule:: Diffusion
#%    
#%  .. |Ito| unicode:: It U+014D
#%

#%  Introduction
#%  ------------
#%  
#%  The functions in this module operate on three conceptual different objects, (although they 
#%  are currently just represented as vectors and arrays.)
#%  
#%  Stochastic processes, denoted ``x``, ``y``, ``w`` are arrays of values which 
#%  are sampled at distance ``dt``, ``ds``, where ``dt``, ``ds`` are either scalar or 
#%  Vectors ``length(dt)=size(W)[end]``.
#%  Stochastic differentials are denoted ``dx``, ``dw`` etc., and are first differences of 
#%  stochastic processes. Finally, t can denote the total time or correspond to a vector of
#%  ``size(W)[end]`` sampling time poins.
#%  
#%  Note the following convention: In analogy with the definition of the |Ito| integral,
#%
#%  	intxdw[i] = x[i]](w[i+1]-w[i]) (== x[i]dw[i])
#%  
#%  and
#%  
#%  	length(w) = length(dw) + 1
#%  




#%  Documentation
#%  -------------
#%  
#%  .. function:: brown1(u, t, n::Integer) 
#%  
#%  	Compute ``n`` equally spaced samples of 1d Brownian motion in
#%  	the interval ``[0,t]``, starting from point ``u``
#%  

function brown1(u, t, n::Integer) 
	x = randn(n)*sqrt(t/(n-1)) ; x[1] = u #(n-1) gaussians
	x = cumsum(x)
end

#%  .. function:: brown1(u, t, n::Integer) 
#%  
#%  	Simulate ``n`` equally spaced samples of ``d``-dimensional Brownian motion in
#%  	the interval ``[0,t]``, starting from point ``u``
#%  


function brown(u, t, d::Integer, n::Integer) 
	x = randn(d, n)*sqrt(t/(n-1)) ; x[:,1] = u #(n-1) gaussians
	x = cumsum(x,2)
end

#%  .. function:: dW1(t, n::Integer)
#%                dW(t, d::Integer, n::Integer)
#%  
#%  	Simulate a ``1``-dimensional (``d``-dimensional)
#%  	Wiener differential with ``n`` values in the 
#%  	the interval ``[0,t]``, starting from point ``u``
#%	


dW1(t, n::Integer) = randn(n)*sqrt(t/n)
dW(t, d::Integer, n::Integer) = randn(d, n)*sqrt(t/n)

#%  .. function:: dW(dt::Vector, d::Integer) 
#%  
#%  	Simulate a ``d``-dimensional Wiener differential sampled at
#%  	time points with distances given by the vector ``dt``		
#%	
mulinto(x,y) = (x[:] = x.*y)
function dW(dt::Vector, d::Integer) 
	#% scale each row of a randn by sqrt(dt)
	dw = Base.each_row(x -> mulinto(x, sqrt(dt)), randn(d, length(dt)))
	return dw
end


### depreciated
dT(t, n::Integer) = ones(n)*(t/n)
function dT0(t, n) 
	dt = ones(n)*(t/(n-1))
	dt[1] = 0
	return dt
end


function dW0(t,n) 
	dw = randn(n)*sqrt(t/(n-1))
	dw[1] = 0
	return dw
end


function dW1(dt) 
	dw = randn(length(dt)).* sqrt(dt)
	return dw
end


#%  .. function:: ito(y, dx)
#%                ito(dx)
#%                cumsum0(dx)
#%  
#%  	Integrate a valued stochastic process with respect to a stochastic differential.
#%  	R, R^2 (d rows, n columns), R^3.
#%  	
#%	``ito(dx)`` is a shortcut for ``ito(ones(size(dx)[end], dx)``.
#%	So ``ito(dx)`` is just a ``cumsum0`` function which is a inverse to ``dx = diff([0, x1, x2, x3,...])``.
#% 
function ito(dx::Vector)
	n = length(dx) + 1
	x = similar(dx, n)
	x[1] = 0
	for i in 2:n
		x[i] += x[i-1] + dx[i-1] 
	end
	x
end

function ito(dx)
	s = [size(dx)...]
	d = ndims(dx)
	n = s[end] += 1
	x = similar(dx, s...)
	if(d==2)
		x[:,1] = 0.0
		for i in 2:n;	x[:,i] += x[:,i-1] + dx[:,i-1]; end	
	elseif (d==3)
		x[:,:,1] = 0.0
		for i in 2:n;	x[:,:,i] += x[:,:,i-1] + dx[:,:,i-1]; end
	else
		error("ndims(dx) > 3 not implemented")
	end
	x
end



sint(dw) = ito(dw)

function ito(x::Matrix, dw::Matrix)
	s = size(dw)[end]
	#determine size ``experimentally''
	dim = [size(x[:,1]*dw[:,1]), n]
	y = similar(dw, dim...)
	y[:,1] = 0
	for i in 2:n
		y[:,1] += y[:,i-1] + x[:,i-1]*dw[:,i-1] 
	end
end


#%  .. function:: ..(y, dx)
#%                ydx(y, dx)
#%  
#%  	``y .. dx'' returns the stochastic differential ``ydx`` defined by the property
#%  
#% 		ito(ydx) == ito(y, dx)
#%  


function ydx(y, dx)
	n = size(dx)[1]
	dim = [n, size(y[1,:]*dx[1,:])...]
	dy = zeros(dim...)
	for i in 2:n+1
		dy[i-1, :] = y[i-1,:]*dx[i-1,:] 
	end
	return dy
end


function ydx(y::Matrix, dx::Matrix)
	n = size(dx)[end]
	dim = [size(y[:,1]'*dx[:,1])..., n]
	dy = zeros(dim...)
	for i in 2:n+1
		dy[:, i-1] = y[:,i-1]'*dx[:, i-1] 
	end
	return dy
end


function (..)(y, dx)
	if (size(y) == ())
		return y*dx
	elseif (size(y)[end] == size(dx)[end])
		return ydx(y,dx)
	else
		return y * dx
	end
end


#%  .. function:: bb(u, v, t, n) 
#%  
#%  	Simulates ``n`` equidistant samples of a Brownian bridge from point ``u`` to ``v`` in time ``t``
#%  

function bb(u, v, t, n) 
	x = cumsum(randn(n)*sqrt(t/(n-1)))
	u + x + dT(v-u-x[end], n)
end

#%  .. function:: dWcond1(v,t,n)
#%  
#%    	Simulates ``n`` equidistant samples of a "bridge noise": that is a Wiener differential ``dW``
#%    	conditioned on ``W(t) = v``
#%  
#%  

function dWcond(v,t,n)
	x = randn(n)*sqrt(t/n)
	x += -sum(x)/n + v/n
end
# sum(dWcond(2, 5, 10)) approx 2

#%  .. function:: aug(dw, dt, n)
#%                aug(dt, n)
#%  
#%  	Take Wiener differential sampled at ``dt`` and return Wiener differential subsampled ``n`` times
#%  	between each observation with new length ``length(dw)*n``.
#%  	``aug(dt,n)`` computes the corresponding subsample of times.
#%  
function aug(dw, dt, n)
	x = zeros(length(dw)*n)
	for i in 1:length(dw)
		x[(1:n)+(i-1)*n] = dWcond(dw[i], dt[i], n)
	end
	x
end

function aug(dt, n)
	x = zeros(length(dt)*n)
	for i in 1:length(dt)
		x[(1:n)+(i-1)*n] = ones(n)*dt[i]/n
	end
	x
end

#todo: aug(dw, dt, ds)

#%  .. function:: quvar(x)
#%               
#%  	Computes quadratic variation of ``x``.
#%  	

function quvar(x)
	sum(diff(x,2).^2)
end

#%  .. function:: bracket(x)
#%                bracket(x,y)
#%  
#%  	Computes quadratic variation process of ``x`` (of ``x`` and ``y``).
#%  	

function bracket(x)
	cumsum0(diff(x).^2)
end

function bracket(x,y)
	cumsum0(diff(x).*diff(y))
end

#%  .. function:: euler(t0, u, b, sigma, dt, dw)
#%                euler(t0, u, b, sigma, dt)
#%  
#%  	Simulates the Euler approximation of a diffusion process
#%  	with drift ``b(t,x)`` and diffusion coefficient ``sigma(t,x)``
#%  	starting in ``(t0, u)`` using ``dt` and given Wiener differential ``dw``.
#%  

#euler: computes the euler approximation of a 
function euler(t0, u, b, sigma, dt, dw)
	X = zeros(length(dw)+1)
	X[1] = u
	t = t0
	for i in 1:length(dw)
		t += dt[i]
		X[i+1] = X[i] + b(t,X[i])*dt[i] + sigma(t,dw[i])*dw[i]
	end
	X
end
euler(t0, u, b, sigma, dt) = euler(t0, u, b, sigma, dt, dW1(dt))
	
end
