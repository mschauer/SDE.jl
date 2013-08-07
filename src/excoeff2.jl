#coefficients for two dimensional example

####### dimension

d = 2

######## starting point

u = [1.,0.]


####### drift 

bB = [ -0.2   -1;   0.5  -0.4]
B = copy(bB)
#beta = [0.5,-0.5]
bbeta = [0., 0.]
beta = copy(bbeta)


b(t,x) = exp(-0.2*t)*bB*x + bbeta
#integral of exp(-0.2*t) to be used 
ph(t,s) = 5.*exp(-0.2*s)-5.*exp(-0.2*t)


####### diffusion coefficient


# uniformly elliptic
function sigma0(s,x) 
 	rho = 2atan(2*x[2])/pi
 	return( 0.7* si* [(2-rho)*[1, 1]  (2+rho)*[1,-1]])
	
end
a0(s,x) = sigma0(s,x)*sigma0(s,x)'

# [0.5,1]-banded elliptic
function sigma1(s,x)
 	m = norm(x-[-0.8, 0.25])
	[1 1; -1  (1 + 1atan(2m)/pi)]
end
a1(s,x) = sigma1(s,x)*sigma1(s,x)'

# only depending on time
sigma2(t,x) = sigma0(0., [sin(t), cos(t)])
a2(s,x) = sigma2(s,x)*sigma2(s,x)'


# choose one
sigma(t,x) = sigma0(t,x)
a(t,x) = a0(t,x)



# fucntion to ignorantly test, if a is uniformly elliptic
function testa(a, N)
	Ra(a,x) = scalar(x'a(0,x)x/(x'x))
	xmax = [NaN, NaN]
	xmin = [NaN, NaN]
	rmax = -Inf
	rmin = Inf
	for i in 1:N
		x = 5*randn(2)
		r = Ra(a, x)
		if (r > rmax)
			rmax = r
			xmax = x
		end 
		if (r < rmin)
			rmin = r
			xmin = x
		end 
		
	end
	return rmin, xmin,rmax, xmax
end
