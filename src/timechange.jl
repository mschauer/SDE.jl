#%  .. function:: tau(s, T)
#%                taui(t, T)
#%  
#%  	Time change mapping s in [0, T] to [0, T] and inverse.
#%  	

tau(s, T) = s.*(2.-s/T) #T*(1.-(1. - s/T)^2)
taui(t, T) =  T-sqrt(T*(T-t))

# h =  T - T*(1.-(1. - s/T)^2) = T(1. - s/T)^2 = T(T - s)^2/T^2 =  (T - s)^2/T 

#%  .. function:: Vtau (s, T, v, B, beta)
#%                dotVtau (s, T, v, B, beta)
#%  
#%  	Time changed V and time changed time derivative of V for generation of U
#%  	


function Vtau (s, T, v, B, beta)
	if (norm(B) < eps2) 
		return v - (T - s)^2/T*beta
	end
 	Binvbeta = B\beta
 	expm(-B*(T - s)^2/T)*( v + Binvbeta) -  Binvbeta
end
function dotVtau (s, T, v, B, beta)
	if (norm(B) < 1E-8) 
		return beta  
	end
 	expm(-B*(T - s)^2/T)*( B*v + beta) 
end

#%  .. function:: UofX(s,x,  T, v,  B, beta)
#%                XofU(s,u,  T, v,  B, beta)
#%                XofU2(S,U, tmin, T, v,  B, beta) 
#%  
#%  Scaled and time changed process U. 
#%  	U(s)= exp(s/2.)*(v(s) - X(tau(s))) 
#%  XofU2 transforms entire process U sampled at time points S.
#%  
	

#careful here, s is in U-time
UofX(s,x,  T, v,  B, beta) = (Vtau(s, T, v, B, beta)- x)/(T-s)

XofU(s,u,  T, v,  B, beta) = Vtau(s, T, v, B, beta)- (T-s)*u

function XofU2(S,U, tmin, T, v,  B, beta) 
 	Ts = similar(S)
	X = similar(U)
	
	for i in 1:length(S)
		s = S[i]
		u = U[:, i]
		Ts[i] = tmin + tau(s,T)
		X[:, i] = Vtau(s, T, v, B, beta)- (T-s)*u
	end
	(Ts,X)

end


#helper functions

ddd(s,u,  T, v,  B, beta) = (tau(s,T), XofU(s,u,  T, v,  B, beta))
ddd(s,u, tmin, T, v,  B, beta) = (tmin+tau(s,T), XofU(s,u,  T, v,  B, beta))


qu(A) = A*A'


	

function J1(t,T, B, A, lambda)
# 	J is responsible for check whether B is zero
	phim = expm(-(T-t)^2/T*B)
	sl = lambda*T/(T-t)^2
	inv( phim*sl*phim'-sl)
end


function J(s,T, B, A, lambda)
	if norm(B) <= eps2
		return lambda
	else
		return J1(s,T, B, A, lambda) 
	end
end


function bU (T, v, b, a, B, beta, lambda)
 (s,x) -> 2./T*dotVtau(s,T,v,B,beta) - 2/T*b(ddd(s, x, T, v, B, beta)...) +  1./(T-s)*(x-   2.*a(ddd(s, x, T, v, B, beta)...)*J(s, T, B, a(T,v), lambda)*x )
end
function bU (tmin, T, v, b, a, B, beta, lambda)
 (s,x) -> 2./T*dotVtau(s,T,v,B,beta) - 2/T*b(ddd(s, x, tmin, T, v, B, beta)...) +   1./(T-s)*(x-   2.*a(ddd(s, x, tmin, T, v, B, beta)...)*J(s, T, B, a(T,v), lambda)*x )
end


function llikeliU(S, U, T, v, b, a,  B, beta, lambda)
	llikeliU(S, U, 0.0, T, v, b, a,  B, beta, lambda)
end

function llikeliU(S, U, tmin, tmax, v, b, a,  B, beta, lambda)
	N = size(U,2)
	T = tmax - tmin
 	function L(s,u)
		j = J(s,T, B, a(tmin+T,v), lambda)
		ju = j*u
		z1 = 2*scalar((b(ddd(s, u, tmin, T, v, B, beta)...)  - B*XofU(s,u,  T, v,  B, beta) - beta)' *ju)
		z2 = -1./(T-s)*trace((a(ddd(s,u, tmin, T, v, B, beta)...) - a(tmin+T,v)) *( j - T*ju*ju' ))
		return  z1 + z2

		
	end
 	som = 0. 
 	for i in 1:N-1
   	  som += L(S[i], leading(U, i))*(S[i+1]-S[i])
   	end
   	som
 	
end


#slightly optimized only computing the likelihood
#replaces consecutive call of euler scheme with drift bu and llikeliU
function transd(S,tmin, tmax, u, v, b, sigma, a, B, beta, lambda, WN, K=1)
#	tmin = T[1]
#	tmax = T[end]
	Delta =tmax - tmin	
	if (norm(B) < eps2) 
		error("case 'norm(B) < eps2' not yet implemented") 
	end
 	Binvbeta = B\beta
 	d = length(u)
 	U = zeros(d, K)
	U[:,1] = exp(0./2)*(expm(-B*Delta*exp(-0.))*( v + Binvbeta) -  Binvbeta - u)

	y = zeros(d)
	z1 = zeros(d)
	jU = zeros(d)
	tBy = zeros(d)
	bUU  = zeros(d)

	for k = 1:K
		U[:,k] = U[:,1]
	end

	Σ = zeros(K)
	N = length(S)
	for i in 1:N-1
		ds = S[i+1]-S[i]
		s = S[i]
		t = Delta*(1.-exp(-s))

		phim = expm(-B*Delta*exp(-s))
		z1[:] = phim *( v + Binvbeta) -  Binvbeta
		j = J(s,Delta, B, a(tmax, v), lambda)
		
		for k in 1:K
			jU[:] = j*U[:,k]
			y[:] = z1  - exp(-s/2)*U[:,k]
			tBy[:] = B*y + beta
			L = (scalar((b(t, y)  - tBy)' *exp(-s/2.)*jU)
			 - 0.5 *trace((a(t, y) - a(tmax,v)) *( j - 1./Delta*(jU)*(jU)' )))
			 #better mapreduce(Multiply(), Add(), a(t, y) - a(tmax,v), j)
			Σ[k] += ds*L
			bUU[:] = Delta*exp(-s/2.)*phim*(B*v + beta)  - Delta*exp(-s/2)*b(t,y) + (0.5*U[:,k] -  a(t, y)*jU)
			U[:, k] = U[:,k] .+  bUU*ds .+ sqrt(Delta)*sigma(t, y)*WN[:,(k-1)*(N-1)+i]*sqrt(ds)
		end	

		#if (i == 10) println("L$i $L ") end

	end
	K > 1 ? Σ : scalar(Σ)
	
end

