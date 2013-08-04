.. currentmodule:: LinProc
  
Introduction to module LinProc
------------------------------

Some functions relevant for generating vector linear processes 

	:math:`dX_t = B X_t + \beta + \sigma d W_t`

where ``B`` is a stable matrix and ``beta`` a vector, ``A = sigma sigma'`` 
and conditional vector linear processes (Ornstein--Uhlenbeck bridges so to say)
ending at time ``T`` in point ``v``,

	:math:`dX^\star_t = B X_t + \beta + A r(s, X^\star_t)  + \sigma d W_t`

where  :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` and 
``p`` is the transition density of ``X``.

The parameter ``lambda`` is the solution to the Lyapunov equation ``B lambda + lambda B' = -A``, see module ``Lyap``, 

     ``lambda = lyap(B', -A)``

If ``B = 0`` set ``lambda = inv(a)``.



Reference 
---------

.. function:: mu(h, x, B, beta)
             
	Expectation :math:`E_x(X_{t})`
	
.. function:: K(h, B, lambda)
             
	Covariance matrix :math:`Cov(X_{t}, X_{t + h})`
	
.. function:: r(h, x, v, B, beta, lambda)
             
	Returns :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` where
	``p`` is the transition density of the linear process, used by ``Bstar``.

.. function:: H(h, B, lambda)
             
	Negative Hessian of :math:`\log p(t,x; T, v) as a function of x.
	
.. function:: bstar(T, v, b, beta, a, lambda)
             
	Returns the drift function of a vector linear process bridge which end at time T in point v.
	
.. function:: bcirc(T, v, b, beta, a, lambda)
             
	Drift for guided proposal derived from a vector linear process bridge which end at time T in point v.
	
.. function:: llikeliXcirc(t, T, Xcirc, b, a,  B, beta, lambda)
             
	Loglikelihood (log weights) of Xcirc with respect to Xstar.
		t, T -- timespan
		Xcirc -- bridge proposal (drift Bcirc and diffusion coefficient sigma) 
		b, sigma -- diffusion coefficient sigma target
		B, beta -- drift b(x) = Bx + beta of Xtilde
		lambda -- solution of the lyapunov equation for Xtilde
.. function:: lp(h, x, y, b, beta, lambda)
             
	Returns :math:`log p(t,x; T, y)`, the log transition density of the linear process, h = T - t 

.. function:: sample_p(h, x, b, beta, lambda) 
             
	Samples from the transition density of the linear process, h = T - t. 

.. function:: lp0(h, x, y,  mu, gamma)
             
	Returns :math:`log p(t,x; T, y)`, the log transition density of a Brownian motion with drift mu and diffusion a=inv(gamma), h = T - t 

.. function:: sample_p0(h, x, mu, l) 
             
	Samples from the transition density a affine Brownian motion. Takes the Cholesky
	factor as argument. 
		l = chol(a)

.. function:: eulerv(t0, u, v, b(s,x), sigma(s,x), Dt, DW::Matrix)
              eulerv(t0, u, b, sigma, dt, dw::Matrix) = eulerv(t0, u, NaN, b, sigma, dt, dw::Matrix)

	Multivariate euler scheme, starting in u, fixing X[N] = v if v!=NaN (this makes sense, if 
	b pulls X towards v.). 	``dw`` -- Wiener differential with ``n`` values in the 
	the interval ``[t0,sum(dt)]`` sampled at timepoints ``t0+Dt[1], t0 + Dt[1] + Dt[2], ...``
	``b, sigma`` -- drift and diffusion coefficient.
	
	Example: 
		
		Dt = diff(linspace(0., T, N))
		DW = randn(2, N-1) .* sqrt(dt)
		dt = Dt[1] yy = euler(0.0, u, b, sigma, Dt, DW)
