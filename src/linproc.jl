# homogeneous vector linear processes with additive noise
module LinProc
#using Randm
export H, r, p, Bstar, Bcirc, Bsharp


#%  .. currentmodule:: LinProc
#%    

#%  Introduction to module LinProc
#%  ------------------------------
#%  
#%  Some functions relevant for generating vector linear processes 
#%  
#%  	:math:`dX_t = B X_t + \beta + \sigma d W_t`
#%  
#%  where ``B`` is a stable matrix and ``beta`` a vector, 
#%  and conditional vector linear processes (Ornstein--Uhlenbeck bridges so to say)
#%  ending at time ``T`` in point ``v``,
#%  
#%  	:math:`dX^\star_t = B X_t + \beta + a r(s, X^\star_t)  + \sigma d W_t`
#%  
#%  where  :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` and 
#%  ``p`` is the transition density of ``X``.
#%  
#%  The parameter ``lambda`` is the solution to the Lyapunov equation ``B lambda + lambda B' = -a``, see module ``Lyap``, 
#%  
#%       ``lambda = lyap(b', -a)``
#%  

#l = lyap(B',-A); B*l + l*B' + A = 0

#%  
#%  Reference 
#%  ---------
#%

#%  .. function:: K(h, lambda, b)
#%               
#%  	Covariance matrix :math:`Cov(X_{t}, X_{t + h})`
#%  	

function B(b, beta)
	(t,x) -> B*x + beta
end	

function SIG(sigma)
	(t,x) -> sigma
end	


function K(h, b, lambda)
	phi = expm(h*b)
	lambda - phi*lambda*phi'
end

function H(h, b, lambda)
	phim = expm(-h*b)
	inv(lambda - phim*lambda*phim')
end

function V(h, v, b, beta)
	binv = inv(b)
	phim = expm(-h*b)
	phim*v + phim * binv * beta - binv * beta 
end
#	binv = inv(B)
#	phim = expm(-h*B)
#	vt =  phim*( v) + phim * binv * beta -binv*beta 


#%  .. function:: r(h, x, v, b, beta, lambda)
#%               
#%  	Returns :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` where
#%  	``p`` is the transition density of the linear process, used by ``Bstar``.
#%  
  

function r(h, x, v, b, beta, lambda)
	H(h, b, lambda)*(x - V(h, v, b, beta))
end

#%  .. function:: Bstar(T, v, b, beta, a, lambda)
#%               
#%  	Returns the drift function of a vector linear process bridge which end at time T in point v.
#%  	

function Bstar(T, v, B, beta, a, lambda)
	(t,x) -> B*x + beta + a * H(T-t, B, lambda)*(x - V(T-t, v, b, beta))
end	

#%  .. function:: Bcirc(T, v, b, beta, a, lambda)
#%               
#%  	Drift for guided proposal derived from a vector linear process bridge which end at time T in point v.
#%  	

function Bcirc(T, v, b, sigma, B, beta, lambda)
	(t,x) -> b(t,x) +  sigma(t,x)*(sigma(t,x))' * H(T-t, B, lambda)*(x - V(T-t, v, B, beta))
end	



function Bsharp(T, v, b )
	(t,x) -> b(t,x) + (v-x)/(T-t)
end


#%  .. function:: lp(h, x, y, b, beta, lambda)
#%               
#%  	Returns :math:`log p(t,x; T, y)`, the log transition density of the linear process, h = T - t 
#%  
function lp(h, x, y, b, beta, lambda)
	z = (x - V(h, y, b, beta))
	(-1/2*length(x)*log(2pi) - 0.5*log(det(K(h,b, lambda))) + 0.5*z'*H(h, b, lambda)*z) 
end
#%  .. function:: p0(h, x, y,  mu, gamma)
#%               
#%  	Returns :math:`log p(t,x; T, y)`, the transition density of a Brownian motion with drift mu and diffusion a=inv(gamma), h = T - t 
#%  
function lp0(h, x, y, mu, gamma)
          (-1/2*length(x)*log(2pi*h) + 0.5*log(det(gamma))  -0.5*(y-x-h*mu)'*gamma*(y-x-h*mu)/h)
         
end

function sample_p0(h, x, mu, l) #l = chol(gamma)
	z = randn(length(x))
	x + l*z*sqrt(h) + h*mu
end



end #linproc
