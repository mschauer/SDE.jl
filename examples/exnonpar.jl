module ExNonpar
using SDE
using SDE.Schauder
require(Pkg.dir("SDE","src", "SDEPlot.jl"))
using SDEPlot


bex1(s,x) = 4(sin(2pi*x) - x)
bex2(s,x) = -4min(x+.5, 0) - 4max(x-.5, 0)
si(s,x) = 1.
range(a) = (minimum(a), maximum(a))


function example(L, logN,T, f, B; VERBOSE = 0)


	global y, t
	srand(42)
	beta = 0.5
	a = -2
	b = 2
	K = 0
	if (B=="B1"); K = 1;end # Schauder basis with constant
	if (B=="B2"); K = 2;end # Schauder basis plus two affine component
				# otherwise Schauder basis

	VERBOSE > 0 && println("L $L, beta $beta")
	N = 10^logN
	dt = T/N
	P = Diffusion{1}(f, si, ())
    tt = linspace(0., T, N+1)
    W = sample(tt, Wiener())
    u = 0.5(b+a)
    y = euler(u, W, P).yy
	
	if (B == "B1")
		truec = fe_transfB1(x-> f(0., x), a,b, L)
	elseif (B == "B2")
		truec = fe_transfB2(x-> f(0., x), a,b, L)
	else
		truec = fe_transf(x-> f(0., x), a,b, L)

	end
    VERBOSE > 0 && println("true drift f.e. coefficients",truec')

	println("range y $(range(y))")

	post = bayes_drift(y, dt, a, b, L, 10ones(K), beta , B) 
 

	pl = visualize_posterior(post,x -> f(0.,x), 1.96)
 
	(post, pl)
end
println("Try: post = ExNonpar.example(5,5, 100, ExNonpar.bex1, \"B2\"); Winston.display(post[2])")

end
