#%  .. function:: tau(s, T)
#%                taui(t, T)
#%  
#%  	Time change mapping s in [0, infty) to [0, T) and inverse.
#%  	

tau(s, T) = T*(1.-exp(-s))
taui(t, T) = -log(1-t/T) 


#%  .. function:: Vtau (s, T, v, B, beta)
#%                dotVtau (s, T, v, B, beta)
#%  
#%  	Time changed V and time changed time derivative of V for generation of U
#%  	

function Vtau (s, T, v, B, beta)
	if (norm(B) < eps2) 
		return v - T*exp(-s)*beta
	end
 	Binvbeta = B\beta
 	expm(-B*T*exp(-s))*( v + Binvbeta) -  Binvbeta
end
function dotVtau (s, T, v, B, beta)
	if (norm(B) < 1E-8) 
		return beta  
	end
 	expm(-B*T*exp(-s))*( B*v + beta) 
end

#%  .. function:: UofX(s,x,  T, v,  B, beta)
#%                XofU(s,u,  T, v,  B, beta)
#%                XofU2(S,U, tmin, T, v,  B, beta) 
#%  
#%  Scaled and time changed process U. 
#%  	U(s)= exp(s/2.)*(v(s) - X(tau(s))) 
#%  XofU2 transforms entire process U sampled at time points S.
#%  
	

UofX(s,x,  T, v,  B, beta) = exp(s/2)*(Vtau(s, T, v, B, beta)- x)

XofU(s,u,  T, v,  B, beta) = Vtau(s, T, v, B, beta)- exp(-s/2)*u

function XofU2(S,U, tmin, T, v,  B, beta) 
 	Ts = similar(S)
	X = similar(U)
	
	for i in 1:length(S)
		s = S[i]
		u = U[1:2, i]
		Ts[i] = tmin + tau(s,T)
		X[1:2, i] = Vtau(s, T, v, B, beta)- exp(-s/2)*u
	end
	(T,X)

end


#helper functions

ddd(s,u,  T, v,  B, beta) = (tau(s,T), XofU(s,u,  T, v,  B, beta))
ddd(s,u, tmin, T, v,  B, beta) = (tmin+tau(s,T), XofU(s,u,  T, v,  B, beta))


qu(A) = A*A'


	

function J1(s,T, B, A, lambda)
# 	J is responsible for check whether B is zero
	phim = expm(-T*exp(-s)*B)
	sl = exp(s)*lambda
	T*inv( phim*sl*phim'-sl)
end

function inv2x2!(j)
	s = 1.0/(j[1,1]*j[2,2]-j[2,1]*j[1,2])
	
	z = j[2,2] 
	j[2,2] = s*j[1,1]
	j[1,1] = s*z

	z = j[2,1] 
	j[2,1] = -s*j[1,2]
	j[2,1] = -s*z
	j
end

#d = 2 numerical stable for s to infty, using putzer's formula for the matrix exponential and then rearranging terms
function J2(s,T, B, A, lambda)
# 	J is responsible for check whether B is zero

	dis = abs2(B[1,1] - B[2,2]) + 4.0*B[1,2]*B[2,1]
	if (dis >= 0.)
		return J1(s,T, B, A, lambda)
	end
 	spur = B[1,1] + B[2,2]

	k = 0.5*spur;
	la = 0.5*sqrt(-dis);

	
	t = -T*exp(-s)	#t = -h, i believe. 
	#solution of -log(-t/T) = s
	
	j = (-exp(2*k*t)/t*(-cos(la*t)*sin(la*t)/la + sin(la*t)^2/la^2*k)*A + # coefficient in this line converges to 1 for s to infty
	-exp(2*k*t)/t*(-(sin(la*t))^2  + (- expm1(-2*k*t) -(cos(la*t)*sinc(la*t/pi)* 2*k*t)) + sin(la*t)^2/la^2* k*k)*lambda +
	-exp(2*k*t)/t*(sin(la*t)^2/la^2*(B*lambda*B')) 
	)

	inv2x2!(j)
#	rearrangement of 
#	I = eye(size(lambda)...)
#	T*exp(-s-2*k*t)* inv(-sin(la*t)^2*lambda +	sin(la*t)^2/(la)^2*(B*lambda*B' + a*k + k*k*I) - cos(la*t)*sin(la*t)/la*(a + 2*k*I))
end

function J(s,T, B, A, lambda)
	if norm(B) <= eps2
		return lambda
	elseif size(lambda) == (2,2) 
		return J2(s,T, B, A, lambda) 
	else
		return J1(s,T, B, A, lambda) 
	end
end


function bU (T, v, b, a, B, beta, lambda)
 (s,x) -> T*exp(-s/2.)*dotVtau(s,T,v,B,beta) - T*exp(-s/2)*b(ddd(s, x, T, v, B, beta)...) + (0.5*eye(size(lambda)...)-  a(ddd(s, x, T, v, B, beta)...)*J(s, T, B, a(T,v), lambda))*x 
end

function bU (tmin, T, v, b, a, B, beta, lambda)
 (s,x) -> T*exp(-s/2.)*dotVtau(s,T,v,B,beta) - T*exp(-s/2)*b(ddd(s, x, tmin, T, v, B, beta)...) + (0.5*eye(size(lambda)...)-  a(ddd(s, x, tmin, T, v, B, beta)...)*J(s, T, B, a(tmin+T,v), lambda))*x 
end



function llikeliU(S, U, T, v, b, a,  B, beta, lambda)
	llikeliU(S, U, 0.0, T, v, b, a,  B, beta, lambda)
end

function llikeliU(S, U, tmin, T, v, b, a,  B, beta, lambda)
	N = size(U,2)
	
 	function L(s,u)
		j = J(s,T, B, a(tmin+T,v), lambda)
		z1 = scalar((b(ddd(s, u, tmin, T, v, B, beta)...)  - B*XofU(s,u,  T, v,  B, beta) - beta)' *j*exp(-s/2.)*u)
		z2 = - 0.5 *trace((a(ddd(s,u, tmin, T, v, B, beta)...) - a(tmin+T,v)) *( j - 1./T*(j*u)*(j*u)' ))
		#z4 =   1./T*0.5 * scalar(u'*j*((a(ddd(s,u, T, v, B, beta)...) - a(T,v))) * j*u)
		#if (abs(z4) > 100) println("Warning: $s $u $z4.") end
		return  z1 + z2

		
	end
 	som = 0. 
	Li1= L(S[1], leading(U, 1))
	for i in 1:N-1
	  Li = Li1
  	  Li1 =	L(S[i+1], leading(U, i+1))
	  som += 0.5*(Li+Li1)*(S[i+1]-S[i])
  	 # if (i == 10) println("Li(1)$i", Li, " ", Li1) end

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

