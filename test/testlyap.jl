
# test lyap with random matrices nmax times for each dimension in dims and sum the difference 
# in norms
function testlyap(dims, nmax)
	tic()
	sum = 0
	for d in dims
		print(d, " ")
		for n = 1:nmax
			# random symmetric matrix with some zeros
			x = randn(d, d) .* randbool(d,d)
			a = x' + x

			# random stable/hurwitz matrix
			x = randn(d, d) .* randbool(d,d)
			a2 = x'*x
			drei =	[ i<j ? 1 : 0 for i in 1:d, j in 1:d]
			y = randn(d, d) .* randbool(d,d)
			
			if(!isposdef(a2)); break;end #test for positive definitesness
		
			b = drei.*y - drei'.*y' - a2
#			print(eigen(b))
			l = atxpxa(b, a)
			sum += norm(b'*l + l*b - a)
		end
	end	
	println("\n Total deviation:", sum)
	toc()
end

function testsyl(dims, nmax)
	tic()
	sum = 0
	for d in dims
		print(d, " ")
		for n = 1:nmax
			# random symmetric matrix
			x = randn(d, d) .* randbool(d,d)
			a = x'+x
			if(!isposdef(a)); break;end #test
			# random stable/hurwitz matrix
			x = randn(d, d) .* randbool(d,d)
			a2 = x'*x
			drei =	[ i<j ? 1 : 0 for i in 1:d, j in 1:d]
			y = randn(d, d) .* randbool(d,d)
			
			if(!isposdef(a2)); break;end #test
		
			b = drei.*y - drei'.*y' - a2
			l = syl(b', b, a)
			sum += norm(b'*l + l*b - a)
		end
	end	
	println("\n Total deviation:", sum)
	toc()
end

