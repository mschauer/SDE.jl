.. currentmodule:: SDE
  
SDE
------- 



Miscellaneous
~~~~~~~~~~~~~
.. function:: syl(a, b, c)
             
    Solves the Sylvester equation ``AX + XB = C``, where ``C`` is symmetric and 
    ``A`` and ``-B`` have no common eigenvalues using (inefficient)
  algebraic approach via the Kronecker product, see http://en.wikipedia.org/wiki/Sylvester_equation


Stochastic Processes
~~~~~~~~~~~~~~~~~~~~

.. function:: mu(t, x, T, P)
         
    Expectation :math:`E_(t,x)(X_{T})`
    
.. function:: K(t, T, P)
         
    Covariance matrix :math:`Cov(X_{T}-x_t)`
    
.. function:: r(t, x, T, v, P)
         
    Returns :math:`r(t,x) = \operatorname{grad}_x \log p(t,x; T, v)` where
    ``p`` is the transition density of the process ``P``.

.. function:: H(t, T, P)
         
    Negative Hessian of :math:`\log p(t,x; T, v)` as a function of ``x``.
    
.. function:: bstar(t, x, T, v, P::MvPro)
         
    Returns the drift function of a vector linear process bridge which end at time T in point v.
    
.. function:: bcirc(t, x, T, v, Pt::Union(MvLinPro, MvAffPro), P::MvPro)
         
    Drift for guided proposal derived from a vector linear process bridge which end at time T in point v.
    
.. function:: lp(t, x, T, y, P)
         
    Returns :math:`log p(t,x; T, y)`, the log transition density of the process ``P``

.. function:: samplep(t, x, T, P) 
         
    Samples from the transition density of the process ``P``.

.. function:: exact(u, tt, P)
         
    Simulate process ``P`` starting in `u` on a discrete grid `tt` from its transition probability.

.. function:: ll(X,P)
         
    Compute log likelihood evaluated in `B`, `beta` and Lyapunov matrix `lambda`
    for a observed linear process on a discrete grid `dt` from its transition density.

.. function:: lp0(h, x, y,  mu, gamma)
         
    Returns :math:`log p(t,x; T, y)`, the log transition density of a Brownian motion with drift mu and diffusion a=inv(gamma), h = T - t 

.. function:: samplep0(h, x, mu, l) 
         
    Samples from the transition density a affine Brownian motion. Takes the Cholesky
    factor as argument. 
        l = chol(a)

.. function:: euler(u, W::CTPath, P::CTPro)

    Multivariate euler scheme for ``U``, starting in ``u`` using the same time grid as the underlying Wiener process ``W``.
    
.. function:: llikeliXcirc(t, T, Xcirc, b, a,  B, beta, lambda)
         
    Loglikelihood (log weights) of Xcirc with respect to Xstar.

        t, T -- timespan
        Xcirc -- bridge proposal (drift Bcirc and diffusion coefficient sigma) 
        b, sigma -- diffusion coefficient sigma target
        B, beta -- drift b(x) = Bx + beta of Xtilde
        lambda -- solution of the lyapunov equation for Xtilde

.. function:: tofs(s, T)
      soft(t, T)

    Time change mapping s in [0, T] (U-time) to t in [t_1, t_2] (X-time), and inverse.
    
.. function:: XofU(UU, tmin, T, v, P) 
  
U is the scaled and time changed process 
    U(s)= exp(s/2.)*(v(s) - X(tofs(s))) 
XofU transforms entire process U sampled at time points ss to X at tt.

.. function:: Vs (s, T, v, B, beta)
      dotVs (s, T, v, B, beta)

    Time changed V and time changed time derivative of V for generation of U
    
.. function:: stable(Y, d, ep)
         
   Return real stable `d`-dim matrix with real eigenvalues smaller than `-ep` parametrized with a vector of length `d*d`, 


    For maximum likelihood estimation we need to search the maximum over all stable matrices.
  These are matrices with eigenvalues with strictly negative real parts.
  We obtain a dxd stable matrix as difference of a antisymmetric matrix and a positive definite matrix.

