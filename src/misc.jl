#miscellenious helper functions and constants

range(x) = (min(x), max(x))


#95% normal quantile

Q95 = sqrt(2)*erfinv(0.95)


eps2 = sqrt(eps())


function scalar(x)
assert(length(x) == 1)
x[1]
end

#compute normal confidence interval for Monte Carlo estimate, precision aware rounding

function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, Q95*ste], int(2-log(10,Q95*ste)))
end

#for sequential estimates, it suffices to provide the running sum Z, running sum of squares Z2 and the number of observations k


function mc2(k, Z, Z2)
	m = Z/k
	
	try
		ste = sqrt(Z2/k-m.^2)/sqrt(k)
	
		res = [m, Q95*ste]
		for i in 1:length(m)
			res[2i-1] = round(res[2i-1], 2- int(log(10,Q95*ste[i])))
			res[2i] = round(res[2i], int(round(2-log(10,Q95*ste[i]))))
		end
		return res
	catch
		return [m, NaN.*m]
	end
	
end
function mc3(k, Z, Z2)
	m = Z/k
	
	try
		va  = sqrt(norm(Z2)/k-norm(m).^2)
		ste = sqrt(Z2/k-m.^2)/sqrt(k)
	
		res = [m, Q95*ste,round(va, 2- int(log(10,va)))]
		for i in 1:length(m)
			res[2i-1] = round(res[2i-1], 2- int(log(10,Q95*ste[i])))
			res[2i] = round(res[2i], int(round(2-log(10,Q95*ste[i]))))
		end
		return res
	catch
		return [m, NaN.*m, NaN]
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
