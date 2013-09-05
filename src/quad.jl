using Cubature
function Hquad(t, T, B, a::Function)
	d = size(B,1)
	function f(s, y)
		y[:] = vec(expm(-(s-t)*B)*a(s)*expm(-(s-t)*B)')
	end
	Q = reshape(hquadrature(d*d, f, t, T; reltol=1E-15, abstol=1E-15, maxevals=5000)[1], d,d) 
	inv(Q)
end	

#obtaining r via quadrature
function varmu(h, x, B, beta)
	function f(s, y)
		y[:] = expm(-s*B)*beta
	end
	integral = Cubature.hquadrature(length(beta), f, 0, h; reltol=1E-15, abstol=1E-15, maxevals=0800)[1]
	expm(h*B)*(x + integral)
end	


function varr(h, x, v, B, beta, lambda)
	mu = varmu(h, x, B, beta)
	expm(h*B')*inv(LinProc.K(h, B, lambda))*(v - mu) 
end

function Qquad(s, T, ph, B, a::Function)
	d = size(B,1)
	
	function f(tau, y)
		y[:] = vec(expm(ph(s,tau)*B)*a(tau)*expm(ph(s,tau)*B)')
	end
	Q = reshape(hquadrature(d*d, f, s, T; reltol=1E-15, abstol=1E-15, maxevals=05000)[1], d,d)
	Q
end
 
function Vquad(s,T, v, ph, B, beta)
	function f(tau, y)
		y[:] = expm(ph(s,tau)*B)*beta(tau)
	end
	expm(-ph(T,s)*B)*v - hquadrature(length(v), f, s, T; reltol=1E-15, abstol=1E-15, maxevals=05000)[1]
	
end


function lp(t, T, x, y, ph, B, beta, a::Function)
 	z = (x -  Vquad(t,T, y, ph, B, t -> beta))
	Q = Qquad(t, T, ph,  B, a)
	l = chol(Q, :L)
	K =  expm(ph(T,t)*B)*Q*expm(ph(T,t)*B)'
	(-1/2*length(x)*log(2pi) -log(apply(*,diag(chol(K)))) - 0.5*norm(l\z)^2) #  - 0.5*log(det(K(h,b, lambda)))
end

