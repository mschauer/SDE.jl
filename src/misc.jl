#miscallenious helper functions

range(x) = (min(x), max(x))

#95% normal quantile

Q95 = sqrt(2)*erfinv(0.95)

#compute normal confidence interval for Monte Carlo estimate, precision aware rounding

function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, Q95*ste], int(2-log(10,Q95*ste)))
end

#for sequential estimates, it suffices to provide the running sum Z, running sum of squares Z2 and the number of observations k

function mc2(K, Z, Z2)
  m = float64(Z/K)
  D = Z2-Z^2/K
  if (D <= eps(m) || K < 2) return ([m, NaN]) end
  ste = float64(sqrt(D)/sqrt((K-1)*K))
  round([m, Q95*ste], int(2-log(10,Q95*ste)))
end


function scalar(x)
	assert(length(x) == 1)
	x[1]
end

function pl(x::Array{ FloatingPoint,2})
	p = FramedPlot()
	#compute range
	R1 = range(x[1,:]) 
	R2 = range(x[2,:])

	#widen a bit	
	R1 = (R1[1] -0.1*(R1[2]-R1[1]),R1[2] + 0.1*(R1[2]-R1[1]))
	R2 = (R2[1]  -0.1*(R2[2]-R2[1]), R2[2] + 0.1*(R2[2]-R2[1]))

	setattr(p, "xrange", R1)
	setattr(p, "yrange", R2)

	add(p, Curve(x[1,:],x[2,:]))
	Winston.display(p)

	p	
end

