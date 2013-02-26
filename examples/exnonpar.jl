require("Itostat")
using Schauder
using Diffusion
using SdeNonparBayes
using Winston


fex1(x) = 4(sin(2pi*x) - x)
fex2(x) = -4min(x+.5, 0) - 4max(x-.5, 0)
range(a) = (min(a), max(a))

function ex1(L, logN,T, f, B)
	global y, t
	srand(42)
	beta = 0.5
	a = -2
	b = 2
	K = 0
	if (B=="B1"); K = 1;end
	if (B=="B2"); K = 2;end

	println("L $L, beta $beta")
	N = 10^logN
	dt = T/N
	t = linspace(0.,T,N)
	y = euler(0, 0.5(b+a), (t,x)-> f(x), (t,x) -> 1, T/N, dW1(T, N))

	t = [0,cumsum(t)]
	if (B == "B1")
		truec = fe_transfB1(x-> f(x), a,b, L)
	elseif (B == "B2")
		truec = fe_transfB2(x-> f(x), a,b, L)
	else
		truec = fe_transf(x-> f(x), a,b, L)

	end
	println("drift f.e. coefficients",truec')
	#pickup!(coeff)
	#println("drift finger coefficients",coeff[finger_pm(L)]')

	println("range y $(range(y))")

	post = bayes_drift(y, dt, a, b, L, 10ones(K), beta , B) 
 

	pl = visualize_posterior(post,f, 1.96)
 
	(post, pl)
end
print("Try: post = ex1(5,5, 100, fex1, \"B2\"); post[2]")
