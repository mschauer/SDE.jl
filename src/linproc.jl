# homogeneous vector linear processes with additive noise
module LinProc
#using Randm
export H, r, p, Bstar, Bcirc, Bsharp, eulerv, llikelixcirc
include("leading.jl")
include("misc.jl")
#%  .. currentmodule:: LinProc
#%    

#%  Introduction to module LinProc
#%  ------------------------------
#%  
#%  Some functions relevant for generating vector linear processes 
#%  
#%  	:math:`dX_t = B X_t + \beta + \sigma d W_t`
#%  
#%  where ``B`` is a stable matrix and ``beta`` a vector, ``A = sigma sigma'`` 
#%  and conditional vector linear processes (Ornstein--Uhlenbeck bridges so to say)
#%  ending at time ``T`` in point ``v``,
#%  
#%  	:math:`dX^\star_t = B X_t + \beta + A r(s, X^\star_t)  + \sigma d W_t`
#%  
#%  where  :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` and 
#%  ``p`` is the transition density of ``X``.
#%  
#%  The parameter ``lambda`` is the solution to the Lyapunov equation ``B lambda + lambda B' = -A``, see module ``Lyap``, 
#%  
#%       ``lambda = lyap(B', -A)``
#%  
#%  If ``B = 0`` set ``lambda = inv(a)``.
#%  
#%  


# l = lyap(B',-A); B*l + l*B' + A = 0

#%  
#%  Reference 
#%  ---------
#%


#%  .. function:: mu(h, x, B, beta)
#%               
#%  	Expectation :math:`E_x(X_{t})`
#%  	
function mu(h, x, B, beta)
	if (norm(B) < eps2) 
		return x + h*beta
	end
	
	binvbeta = B\beta
	phi = expm(h*B)
	phi*(x + binvbeta) - binvbeta
end	

#%  .. function:: K(h, B, lambda)
#%               
#%  	Covariance matrix :math:`Cov(X_{t}, X_{t + h})`
#%  	

function K(h, B, lambda)
	if (norm(B)) < eps2 
		return h*inv(lambda)
	end
	phi = expm(h*B)
	lambda - phi*lambda*phi'
end

#%  .. function:: r(h, x, v, B, beta, lambda)
#%               
#%  	Returns :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` where
#%  	``p`` is the transition density of the linear process, used by ``Bstar``.
#%  
  

function r(h, x, v, B, beta, lambda)
	H(h, B, lambda)*(V(h, v, B, beta)-x)
end



#%  .. function:: H(h, B, lambda)
#%               
#%  	Negative Hessian of :math:`\log p(t,x; T, v) as a function of x.
#%  	

function H(h, B, lambda)
	if (norm(B)) < eps2 
		return lambda/h
	end
	phim = expm(-h*B)
	inv(phim*lambda*phim'-lambda)
end


function L(h, B, lambda)
	if (norm(B)) < eps2 
		return chol(h*inv(lambda))
	end
	phim = expm(-h*B)
	chol(phim*lambda*phim'-lambda, :L)
end

# x'inv(K)*x =  norm(chol(K, :L)\x)^2

# technical function

function V(h, v, B, beta)
	if (norm(B) < eps2) 
		return v - h*beta
	end
	binvbeta = B\beta
	phim = expm(-h*B)
	phim*(v + binvbeta) - binvbeta  #CHECK
end


#%  .. function:: bstar(T, v, b, beta, a, lambda)
#%               
#%  	Returns the drift function of a vector linear process bridge which end at time T in point v.
#%  	

function bstar(T, v, B, beta, a, lambda)
	(t,x) -> B*x + beta + a * H(T-t, B, lambda)*(V(T-t, v, b, beta)-x)
end	

#%  .. function:: bcirc(T, v, b, beta, a, lambda)
#%               
#%  	Drift for guided proposal derived from a vector linear process bridge which end at time T in point v.
#%  	

function bcirc(T, v, b, sigma, B, beta, lambda)
	(t,x) -> b(t,x) +  sigma(t,x)*(sigma(t,x))' * H(T-t, B, lambda)*(V(T-t, v, B, beta)-x)
end	


#%  .. function:: llikeliXcirc(t, T, Xcirc, b, a,  B, beta, lambda)
#%               
#%  	Loglikelihood (log weights) of Xcirc with respect to Xstar.
#%  		t, T -- timespan
#%  		Xcirc -- bridge proposal (drift Bcirc and diffusion coefficient sigma) 
#%  		b, sigma -- diffusion coefficient sigma target
#%  		B, beta -- drift b(x) = Bx + beta of Xtilde
#%  		lambda -- solution of the lyapunov equation for Xtilde

function llikeliXcirc(t, T, Xcirc, b, a,  B, beta, lambda)
	N = size(Xcirc,2)
	v = leading(Xcirc, N) #like [X, n]
	function L(s,x)
		R = LinProc.H(T-s, B, lambda)*(LinProc.V(T-s, v, B, beta)-x)
	  	return  (b(s,x) - B*x - beta)' * R - 0.5 *trace((a(s,x) - a(T,v)) *(LinProc.H(T-s, B, lambda) - R*R'))
	end
	
	sum = 0
	s= t
	x = similar(v)
	for i in 0:N-1-1 #skip last value, summing over n-1 elements
	  s = t + (T-t)*(i)/(N-1) 
	  x = leading(Xcirc, i+1)
	  sum += scalar(L(s, x)) * (T-t)/N
	end
	sum += scalar( 2*sqrt(T-s)*(b(T,x) - B*x - beta)'* a(T,v)*(v-x)) #interpolate drift part of last interval like square root
	
	sum
end



#%  .. function:: lp(h, x, y, b, beta, lambda)
#%               
#%  	Returns :math:`log p(t,x; T, y)`, the log transition density of the linear process, h = T - t 
#%  
function lp(h, x, y, B, beta, lambda)
	if (norm(B) < eps2)
		return lp0(h, x, y, beta, lambda)
	end
	z = (x - V(h, y, B, beta))
	l = L(h, B, lambda)
	(-1/2*length(x)*log(2pi) -log(apply(*,diag(chol(K(h,B, lambda))))) - 0.5*norm(l\z)^2) #  - 0.5*log(det(K(h,b, lambda)))
end

#%  .. function:: sample_p(h, x, b, beta, lambda) 
#%               
#%  	Samples from the transition density of the linear process, h = T - t. 
#%  

function sample_p(h, x, B, beta, lambda) 
	if (norm(B) < eps2)
		z = randn(length(x))
		return x + chol(lambda)\z*sqrt(h) + h*beta

	end
	phi = expm(h*B)
	binvbeta = B\beta

	mu = phi*(x + binvbeta) - binvbeta 
	k = lambda - phi*lambda*phi'
	l = chol(k)

	z = randn(length(x))
	mu + l*z
end


#%  .. function:: lp0(h, x, y,  mu, gamma)
#%               
#%  	Returns :math:`log p(t,x; T, y)`, the log transition density of a Brownian motion with drift mu and diffusion a=inv(gamma), h = T - t 
#%  
function lp0(h, x, y, mu, gamma)
          (-1/2*length(x)*log(2pi*h) + 0.5*log(det(gamma))  -0.5*(y-x-h*mu)'*gamma*(y-x-h*mu)/h)[1]
         
end
#%  .. function:: sample_p0(h, x, mu, l) 
#%               
#%  	Samples from the transition density a affine Brownian motion. Takes the Cholesky
#%  	factor as argument. 
#%  		l = chol(a)
#%  

function sample_p0(h, x, mu, l) #l = chol(a)
	z = randn(length(x))
	x + l*z*sqrt(h) + h*mu
end


#%  .. function:: eulerv(t0, u, v, b(s,x), sigma(s,x), Dt, DW::Matrix)
#%                eulerv(t0, u, b, sigma, dt, dw::Matrix) = eulerv(t0, u, NaN, b, sigma, dt, dw::Matrix)
#%  
#%  	Multivariate euler scheme, starting in u, fixing X[N] = v if v!=NaN (this makes sense, if 
#%  	b pulls X towards v.). 	``dw`` -- Wiener differential with ``n`` values in the 
#%  	the interval ``[t0,sum(dt)]`` sampled at timepoints ``t0+Dt[1], t0 + Dt[1] + Dt[2], ...``
#%	``b, sigma`` -- drift and diffusion coefficient.
#%  	
#%  	Example: 
#%  		
#%  		Dt = diff(linspace(0., T, N))
#%  		DW = randn(2, N-1) .* sqrt(dt)
#%  		dt = Dt[1] yy = euler(0.0, u, b, sigma, Dt, DW)
		

eulerv(t0, u, b, sigma, dt, dw::Matrix) = eulerv(t0, u, NaN, b, sigma, dt, dw::Matrix)

function eulerv(t0, u, v, b, sigma, dt, dw::Matrix)
	S = size(dw)
	N = S[end] + 1
	
	shape = size(sigma(0,u)*leading(dw,1))
 
	X = zeros(shape..., N)

	y = copy(u)
	t = t0

	
	for i in 1:N-1
		subleading(X,i)[:] = y
		t += dt[i]
		y[:] = y .+  b(t,y)*(dt[i]) .+ sigma(t,y)*leading(dw, i)
	
	end
	
	if (v == NaN)
		subleading(X,N)[:] = v
	else
		subleading(X,N)[:] = y
	end
	X
end


include("clark.jl")

end #linproc
