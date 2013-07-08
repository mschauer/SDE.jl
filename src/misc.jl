#improvised misc
range(x) = (min(x), max(x))

function scalar(x)
assert(length(x) == 1)
x[1]
end

function mc(M)
  m = mean(M)
  ste = std(M)/sqrt(length(M))
  round([m, 1.96*ste], int(2-log(10,1.96*ste)))
end

function mc2(k, LL, LL2)
	m = LL/k
	
	try
		ste = sqrt(LL2/k-m.^2)/sqrt(k)
	
		res = [m, 1.96*ste]
		for i in 1:length(m)
			res[2i-1] = round(res[2i-1], 2- int(log(10,1.96*ste[i])))
			res[2i] = round(res[2i], int(round(2-log(10,1.96*ste[i]))))
		end
		return res
	catch
		return [m, NaN.*m]
	end
	
end
