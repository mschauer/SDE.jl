using Lyap

# test lyap with random matrices nmax times for each dimension in dims and sum the difference 
# in norms
function testlyap(dims, n)
	tic()
	sum = 0
	for d in dims
		print("dims = $dims, n = $n")
		for k = 1:n
			# random symmetric matrix with some zeros
			x = randn(d, d) .* randbool(d,d)
			a = x' + x

			# random stable/hurwitz matrix
			a2 = a
			while true
			x = randn(d, d) .* randbool(d,d)
			a2 = x'*x
			if(!isposdef(a2))
				print(".")
			else 
				break
			end #test for positive definitesness
			
			
			end
			drei =	[ i<j ? 1 : 0 for i in 1:d, j in 1:d]
			y = randn(d, d) .* randbool(d,d)
			
			
#			b = drei.*y - drei'.*y' - a2
			println(eig(b)[1])
			l = lyap(b, a)
			sum += norm(b'*l + l*b - a)
		end
	end	
	println("\nTotal deviation:", sum)
	toc()
end

# test syl with random matrices n times for each dimension in dims and sum the difference 
# in norms
function testsyl(dims, n)
	tic()
	sum = 0
	for d in dims
		print("dims = $dims, n = $n")
		for k = 1:n
			# random symmetric matrix
			x = randn(d, d) .* randbool(d,d)
			a = x'+x
			# random stable/hurwitz matrix
			a2 = a
			while true
			x = randn(d, d) .* randbool(d,d)
			a2 = x'*x
			if(!isposdef(a2))
				print(".")
			else 
				break
			end #test for positive definitesness
			
			
			end
			
			
			drei =	[ i<j ? 1 : 0 for i in 1:d, j in 1:d]
			y = randn(d, d) .* randbool(d,d)
			
				
			b = drei.*y - drei'.*y' - a2
			l = syl(b', b, a)
			sum += norm(b'*l + l*b - a)
		end
	end	
	println("\nTotal deviation:", sum)
	toc()
end

