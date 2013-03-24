module Lyap
export lyap, syl

#%  .. currentmodule:: Lyap
#%    
#%  Introduction to module Lyap
#%  ----------------------------
#%        
#%   
#%    DESCRIPTION
#%         Solves the real matrix equation A'X + XA = C, where A and C are
#%         constant matrices of dimension n x n with C=C'.  The matrix A
#%         is transformed into upper Schur form and the transformed system
#%         is solved by back substitution.  The option is provided to input
#%         the Schur form directly and bypass the Schur decomposition.
#%         This equation is also know as continous Lyapunov equation.
#%         
#%         The method of Bartels and Stewart is used. 
#%         The system is first reduced such that A is in upper real schur
#%         form. The resulting triangular system is solved via back-substitution.
#%         Has a unique solution, if A and -A have no common eigenvalues, which is guaranteed
#%         if A is stable (and the real part of each eigenvalue is negative).
#%   
#%    HISTORY
#%         The classic ACM algorithm from Bartels and Stewart was implemented E. Armstrong 
#%         as part of ORACLS -- optimal regulator algorithms for the control of linear systems.
#%         The implementation from the nasa cosmic archive is reported to be in the 
#%         public domain, under the terms of Title 17, Chapter 1, 
#%         Section 105 of the US Code. This is rather direct translation of forementioned
#%         implementation into Julia, put under MIT licence as Julia.
#%   
#%    REFERENCES
#%         * Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
#%           the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
#%           Sept. 1972, pp. 820-826.
#%   
#%    SEE ALSO
#%         atxpxa in ORACLS, strsyl, dtrsyl in LAPACK, lyap in GNU Octave
#%  
#%  
#%  

#%  Reference 
#%  ---------
#%  
#%  .. function:: issquare(a)
#%               
#%  	Checks if matrix ``a`` is square.
#%  	


function issquare(a::Matrix)
	 size(a,2) == size(a,1)
end


#%  .. function:: lyap(a, cc)
#%  
#%         Solves the real matrix equation A'X + XA = C, where A and C are
#%         constant matrices of dimension n x n with C=C'.  The matrix A
#%         is transformed into upper Schur form and the transformed system
#%         is solved by back substitution.  The option is provided to input
#%         the Schur form directly and bypass the Schur decomposition.
#%         This equation is also know as continous Lyapunov equation.
#%         
#%         The method of Bartels and Stewart is used. 
#%         The system is first reduced such that A is in upper real schur
#%         form. The resulting triangular system is solved via back-substitution.
#%         Has a unique solution, if A and -A have no common eigenvalues, which is guaranteed
#%         if A is stable (and the real part of each eigenvalue is negative).
#%  

function lyap (a, cc)
	assert(issquare(a))
	assert(issquare(cc))
	c = copy(cc)
	a, u = schur(a) # u*a*u'
	n = size(a,1)
	r = similar(a, n) # A(N1,j)
	s = similar(a, n) # A(i, N1)
	# transform c
	for i in 1:n
		c[i,i] /= 2.
	end
	for i in 1:n
		for j in 1:n
			r[j] = 0.
	        for k in i:n
            	r[j] += c[i,k]*u[k,j]
            end
        end
        for j in 1:n
          c[i,j] = r[j]
		end
	end
	for j in 1:n
		for i in 1:n
		s[i] = 0.
			for k in 1:n
				s[i] += u[k,i]*c[k, j]
			end
		end
		for i in 1:n
			c[i,j] = s[i]
		end    
	end
	for i in 1:n
		for j in i:n
			c[i,j] += c[j,i]
			c[j,i] = c[i,j]
		end
	end

	c = symslv(a, c)
	
	
	# transform c back to the solution.
	for i in 1:n
		c[i,i] /= 2.
	end

	for i in 1:n
		for j in 1:n
          r[j] = 0.
          for k in i:n
            r[j] += c[i,k]*u[j,k]
		  end
		end
		for j in 1:n
			  c[i,j] = r[j]
		end
	end
	for j in 1:n
		for i in 1:n
          s[i] = 0.
          for k in 1:n
            	s[i] += + u[i,k]*c[k,j]
		  end
		end
		for i in 1:n
          c[i,j] = s[i]
	  	end
	end
	for i in 1:n
		for j in i:n
			c[i,j] += c[j,i]
			c[j,i] = c[i,j]
		end
	end
	c
end


#%  .. function:: symslv(a, c) 
#%               
#%      Solves ``A'*x + x*A = C``, where ``C`` is symmetric and ``A`` is in upper real schur form.
#%      via back substitution
#%  

function symslv(a, c) 
	t = similar(a, 4,4)
	t3 = similar(a, 3,3)
	t2 = similar(a, 2,2)
	
	p = similar(a, 5)
	n = size(a,1)
	l = 1
	while true #10
	
		dl = 1
		if l != n
			if a[l+1,l] != 0.
				dl = 2
			end
		end
	
		ll = l+dl-1
		k = l
		while k <= n 
			km1 = k-1
			dk = 1
			if k != n
				if a[k+1,k] != 0. 
					dk = 2
				end
			end
			
			kk = k+dk-1
			if k != l
				for i=k:kk
				for j=l:ll
				for ia=l:km1
					c[i,j] = c[i,j] - a[ia,i]*c[ia,j]
				end
				end
				end
			end
			
			
			if dl != 2 
				if dk != 2 
					tt = a[k,k] + a[l,l]
					if tt == 0.0
						error("division by zero")
					end 
					c[k,l] = c[k,l]/tt #try
				else
					t2[1,1] = a[k,k] + a[l,l]
					t2[1,2] = a[kk,k]
					t2[2,1] = a[k,kk]
					t2[2,2] = a[kk,kk] + a[l,l]
					p = [c[k,l], c[kk,l]]
					p = t2\p
					c[k,l] = p[1]
					c[kk,l] = p[2]
				end
			elseif dk != 2
				t2[1,1] = a[k,k] + a[l,l]
				t2[1,2] = a[ll,l]
				t2[2,1] = a[l,ll]
				t2[2,2] = a[k,k] + a[ll,ll]
				p = [c[k,l], c[k,ll]]
				p = t2\p
				c[k,l] = p[1]
				c[k,ll] = p[2]

			elseif k == l 
				t3[1,1] = a[l,l]
				t3[1,2] = a[ll,l]
				t3[1,3] = 0.
				t3[2,1] = a[l,ll]
				t3[2,2] = a[l,l] + a[ll,ll]
				t3[2,3] = t3[1,2]
				t3[3,1] = 0.
				t3[3,2] = t3[2,1]
				t3[3,3] = a[ll,ll]
				p = [c[l,l]/2,  c[ll,l],  c[ll,ll]/2]
				p = t3\p[1:3]
				c[l,l] = p[1]
				c[ll,l] = p[2]
				c[l,ll] = p[2]
				c[ll,ll] = p[3]
			else
				t[1,1] = a[k,k] + a[l,l]
				t[1,2] = a[kk,k]
				t[1,3] = a[ll,l]
				t[1,4] = 0.
				t[2,1] = a[k,kk]
				t[2,2] = a[kk,kk] + a[l,l]
				t[2,3] = 0.
				t[2,4] = t[1,3]
				t[3,1] = a[l,ll]
				t[3,2] = 0.
				t[3,3] = a[k,k] + a[ll,ll]
				t[3,4] = t[1,2]
				t[4,1] = 0.
				t[4,2] = t[3,1]
				t[4,3] = t[2,1]
				t[4,4] = a[kk,kk] + a[ll,ll]
				p = [c[k,l], c[kk,l], c[k,ll], c[kk,ll]]
				p = t\p
				c[k,l] = p[1]
				c[kk,l] = p[2]
				c[k,ll] = p[3]
				c[kk,ll] = p[4]
			end
			k = k + dk
		end #while k <= n
	
		ldl = l + dl
		if ldl > n
			return c
		end
		
		for j = ldl:n
			for i = l:ll
				c[i,j] = c[j,i]
			end
			for i = j:n
				for k = l:ll
					c[i,j] = c[i,j] - c[i,k]*a[k,j] - a[k,i]*c[k,j]
				end
				c[j,i] = c[i,j]
			end
		end
		l = ldl
	end #10
end




#%  .. function:: syl(a, b, c)
#%               
#%  	Solves the Sylvester equation ``AX + XB = C``, where ``C`` is symmetric and 
#%  	``A`` and ``-B`` have no common eigenvalues using (inefficient)
#%	algebraic approach via the Kronecker product, see http://en.wikipedia.org/wiki/Sylvester_equation
#%  

function syl(a, b, c)
	assert(issquare(a))
	assert(issquare(b))
	assert(issquare(c))

	k = kron(eye(a), a) + kron(b', eye(b))
	xvec=k\vec(c)
	reshape(xvec,size(c))
end




end
