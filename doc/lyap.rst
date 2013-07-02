.. currentmodule:: Lyap
  
Introduction to module Lyap
----------------------------
      
 
  DESCRIPTION
       Solves the real matrix equation A'X + XA = C, where A and C are
       constant matrices of dimension n x n with C=C'.  The matrix A
       is transformed into upper Schur form and the transformed system
       is solved by back substitution.  The option is provided to input
       the Schur form directly and bypass the Schur decomposition.
       This equation is also know as continuous Lyapunov equation.
       
       The method of Bartels and Stewart is used. 
       The system is first reduced such that A is in upper real schur
       form. The resulting triangular system is solved via back-substitution.
       Has a unique solution, if A and -A have no common eigenvalues, which is guaranteed
       if A is stable (and the real part of each eigenvalue is negative).
 
  HISTORY
       The classic ACM algorithm from Bartels and Stewart was implemented E. Armstrong 
       as part of ORACLS -- optimal regulator algorithms for the control of linear systems.
       The implementation from the nasa cosmic archive is reported to be in the 
       public domain, under the terms of Title 17, Chapter 1, 
       Section 105 of the US Code. This is rather direct translation of forementioned
       implementation into Julia, put under MIT licence as Julia.
 
  REFERENCES
       * Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
         the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
         Sept. 1972, pp. 820-826.
 
  SEE ALSO
       atxpxa in ORACLS, strsyl, dtrsyl in LAPACK, lyap in GNU Octave



Reference 
---------

.. function:: issquare(a)
             
	Checks if matrix ``a`` is square.
	
.. function:: lyap(a, cc)

       Solves the real matrix equation A'X + XA = C, where A and C are
       constant matrices of dimension n x n with C=C'.  The matrix A
       is transformed into upper Schur form and the transformed system
       is solved by back substitution.  The option is provided to input
       the Schur form directly and bypass the Schur decomposition.
       This equation is also know as continuous Lyapunov equation.
       
       The method of Bartels and Stewart is used. 
       The system is first reduced such that A is in upper real schur
       form. The resulting triangular system is solved via back-substitution.
       Has a unique solution, if A and -A have no common eigenvalues, which is guaranteed
       if A is stable (and the real part of each eigenvalue is negative).

.. function:: symslv(a, c) 
             
    Solves ``A'*x + x*A = C``, where ``C`` is symmetric and ``A`` is in upper real schur form.
    via back substitution

.. function:: syl(a, b, c)
             
	Solves the Sylvester equation ``AX + XB = C``, where ``C`` is symmetric and 
	``A`` and ``-B`` have no common eigenvalues using (inefficient)
	algebraic approach via the Kronecker product, see http://en.wikipedia.org/wiki/Sylvester_equation

