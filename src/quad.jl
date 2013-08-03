using Cubature
function Hquad(t, T, B, a::Function)
	function f(s, v)
		v[:] = vec(expm(-(s-t)*B)*a(s)*expm(-(s-t)*B)')
	end
	Q = reshape(hquadrature(d*d, f, t, T, 1E-15, 1E-15, 5000)[1], d,d)
	inv(Q)
end	

#obtaining r via quadrature
function varmu(h, x, B, beta)
	function f(s, y)
		y[:] = expm(-s*B)*beta
	end
	integral = Cubature.hquadrature(d, f, 0, h, 1E-15, 1E-15, 800)[1]
	expm(h*B)*(x + integral)
end	


function varr(h, x, v, B, beta, lambda)
	mu = varmu(h, x, B, beta)
	expm(h*B')*inv(LinProc.K(h, B, lambda))*(v - mu) 
end

function Qquad(s, T, ph, B, a::Function)
	
	function f(tau, v)
		v[:] = vec(expm(ph(s,tau)*B)*a(tau)*expm(ph(s,tau)*B)')
	end
	Q = reshape(hquadrature(d*d, f, s, T, 1E-15, 1E-15, 5000)[1], d,d)
	Q
end
 
function Vquad(s,T, v, ph, b, beta)
	function f(tau, y)
		y[:] = expm(ph(s,tau)*b)*beta(tau)
	end
	expm(-ph(T,s)*b)*v - hquadrature(d, f, s, T, 1E-15, 1E-15, 5000)[1]
	
end


function lp(t, T, x, y, ph, B, beta, a::Function)
 	z = (x -  Vquad(t,T, y, ph, B, t -> beta))
	Q = Qquad(t, T, ph,  B, a)
	l = chol(Q, :L)
	K =  expm(ph(T,t)*B)*Q*expm(ph(T,t)*B)'
	(-1/2*length(x)*log(2pi) -log(apply(*,diag(chol(K)))) - 0.5*norm(l\z)^2) #  - 0.5*log(det(K(h,b, lambda)))
end

