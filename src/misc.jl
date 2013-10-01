#miscellenious helper functions and constants

range(x) = (min(x), max(x))

#95% normal quantile

Q95 = sqrt(2)*erfinv(0.95)

roundv(v::Vector, i) = map(i -> round(i, 2), v)
roundv(v , i) = round(v,i)

eps2 = sqrt(eps())


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
  min(max(a, x), b)
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


function isignif_og(x, digits, base)
    if base == 10
        ifloor(log10(abs(x)) - digits + 1.)
    elseif base == 2
        ifloor(log2(abs(x)) - digits + 1.)
    else
        ifloor(log2(abs(x))/log2(base) - digits + 1.)
    end
end

function isignif(x, digits::Integer, base::Integer=10)
    if digits < 0
        throw(DomainError())
    end
    if x==0 || !isfinite(x)
        return x, 0
    end
    og = isignif_og(float(x), digits, base)
    iround(float(x)/float(base)^og), og
end


strisignif(xi, og) = string(xi * og)

function strisignif(xi::Integer, og::Integer)
  
	if og >= 0 && xi == 0
		"0"
	elseif og >= 0 && xi != 0
		string(xi) * "0"^(og)
	elseif abs(xi) >= 10^(-og) 
		s = string(xi)
 		s[1: end+og]*"."*s[end+og+1:end]
	else
		s = string(abs(xi))
		(xi < 0 ? "-":"")* "0."*"0"^(-length(s)-og)*s
	end
end


function roundste(m, ste::FloatingPoint)

 if (m == -0.0) m = 0.0 end
 if (ste == -0.0) ste = 0.0 end

 if 
 	!isfinite(m) return string(m)
 end
 
 if isnan(ste) 
 	return (isinteger(m) ?  string(int(m)) : string(m)) * " ± NaN (se)" 
 elseif isinf(ste)
 	return (isinteger(m)  ? string(int(m)) : string(m)) * " ± Inf (se)" 
 end
 assert(ste >= 0.) 
 stei, og = isignif(ste, 3)
 mi =  iround(float(m)/float(10)^og)
 
 strisignif(mi, og)*" ± "*strisignif(stei, og)*(" (se)")
end

function roundste(r)
 	m, ste = r
 	res = roundste(m[1], ste[1]) 
	for i in 2:length(m)
		res *= ", "* roundste(m[i], ste[i]) 
	end
	res
end


strmc2(k, Z, Z2) = roundste( mc2(k, Z, Z2, false))


#stdv(k, Z, Z2) = va(k, Z, Z2)/k
# X = [ZL, L]
# X2 = X * X'	

function selfadjmc(k, X, X2)
	m = X[1]/X[2]
	v = (va(k, X[1], X2[1,1]) - 2m*(X2[1,2]/k - X[1]*X[2]/k^2) + m^2*va(k, X[2], X2[2,2]))/k
#	println([va(k, X[1], X2[1,1]), - 2m*(X2[1,2]/k - X[1]*X[2]/k^2),m^2*va(k, X[2], X2[2,2])])
	try
		ste = sqrt(v)
		res = [m, ste]
		return res
	catch
		return [m, NaN.*m]
	end
	
end

#plot (2d) sample path

function pl(x::Array{ FloatingPoint,2})
	p = FramedPlot()
	#compute range
	R1 = range(x[1,:]) 
	R2 = range(x[2,:])

	#widen a bit	
	R1 = (R1[1] -0.1*(R1[2]-R1[1]), R1[2] + 0.1*(R1[2]-R1[1]))
	R2 = (R2[1] -0.1*(R2[2]-R2[1]), R2[2] + 0.1*(R2[2]-R2[1]))

	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

	add(p, Curve(x[1,:],x[2,:]))
	Winston.display(p)

	p	
end
