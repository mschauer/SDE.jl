.. currentmodule:: LinProc
  
Introduction to module LinProc
------------------------------

Some functions relevant for generating vector linear processes 

	:math:`dX_t = B X_t + \beta + \sigma d W_t`

where ``b`` is a stable matrix and ``beta`` a vector, 
and conditional vector linear processes (Ornstein--Uhlenbeck bridges so to say)
ending at time ``T`` in point ``v``,

	:math:`dX^\star_t = B X_t + \beta + a r(s, X^\star_t)  + \sigma d W_t`

where  :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` and 
``p`` is the transition density of ``X``.

The parameter ``lambda`` is the solution to the Lyapunov equation ``B'lambda + lambda B = -a``, see module ``Lyap``, 

     ``lambda = lyap(b', -a)``


Reference 
---------

.. function:: K(h, lambda, b)
             
	Covariance matrix :math:`Cov(X_{t}, X_{t + h})`
	
.. function:: r(h, x, v, b, beta, lambda)
             
	Returns :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` where
	``p`` is the transition density of the linear process, used by ``Bstar``.

.. function:: Bstar(T, v, b, beta, a, lambda)
             
	Returns the drift function of a vector linear process bridge which end at time T in point v.
	
.. function:: Bstar(T, v, b, beta, a, lambda)
             
	Drift for guided proposal derived from a vector linear process bridge which end at time T in point v.
	
