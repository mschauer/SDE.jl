 
d = 2

bB = [ -0.2   -1;   0.5  -0.4]
B = copy([ -0.2   -1;   0.5  -0.4])
#beta = [0.5,-0.5]
beta = [0., 0.]

# uniformly elliptic
function sigma0(s,x) 
	m = norm(x)
	if (m > eps())  
		rho = 1 + 5*2atan(m)/pi
		return( si*((1-2atan(m)/pi)*eye(2)+ 2atan(m)/pi/m*[[x[2], -x[1]]  [rho*x[1], rho*x[2]]]))
	end
	return si*eye(2)
end

# banded elliptic
function sigma1(s,x)
 	m = norm(x-u)
	[1 1; -1  (1 + 8atan(2m)/pi)]
end


#sigma(t,x) = sigma0(0., [sin(t), cos(t)])
#sigma(t,x) = sigma1(t,x)
sigma(t,x) = sigma0(t,x)
a(s,x) = sigma(s,x)*sigma(s,x)'



b(t,x) = exp(-0.2*t)*bB*x + beta
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)
u = [1.,0.]
