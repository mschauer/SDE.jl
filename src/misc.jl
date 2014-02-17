#miscellaneous helper functions and constants
Q95 = sqrt(2.)*erfinv(0.95)


range(x) = (minimum(x), maximum(x))

#95% normal quantile


roundv(v::Vector, i) = map(i -> round(i, 2), v)
roundv(v , i) = round(v,i)

eps2 = sqrt(eps())




intervalgaussian(U, a, b) = broadcast(intervalgaussian, U, a, b)
function intervalgaussian (U::Real,a::Real, b::Real)
	a, b = min(a,b), max(a,b)
	c =  erf(a/sqrt(2.))   
	d =  erf(b/sqrt(2.))  
	sqrt(2.)*erfinv(c + U.*(d-c))  
end



function scalar(x)
assert(length(x) == 1)
x[1]
end

function norma!(x)
 minx = min(x)
 x = (x - minx) / (max(x)-minx)
 x
end


function cut(x, a, b)
	warn("cut(): better use clamp()")
	clamp(x, a,b)
end

#compute normal confidence interval for Monte Carlo estimate, precision aware rounding

function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, Q95*ste], int(2-log(10,Q95*ste)))
end

#for sequential estimates, it suffices to provide the running sum Z, running sum of squares Z2 and the number of observations k


va(k, Z, Z2) =  Z2/k-(Z/k)^2
ste(k, Z, Z2) = sqrt(Z2/k^2-(Z/k)^2/k)



function mc2(k, Z, Z2, normalci = true)
	m = Z/k

	if normalci F = Q95 else F = 1. end

	try
		ste = sqrt(Z2/k-m.^2)/sqrt(k)
	
		res = (m, F*ste)
	
		
		return res
	catch
		return (m, NaN.*m)
	end
	
end


function selfadjmc(k, X, X2)
	m = X[1]/X[2]
	v = (va(k, X[1], X2[1,1]) - 2m*(X2[1,2]/k - X[1]*X[2]/k^2) + m^2*va(k, X[2], X2[2,2]))/k
	try
		ste = sqrt(v)
		res = [m, ste]
		return res
	catch
		return [m, NaN.*m]
	end
	
end

