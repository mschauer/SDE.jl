include("StoPro.jl")
using Distributions
using StoPro
srand(8)
# sum([ norm(mean([ sample(Wiener(d),  linspace(0.,2.,5)).X[1,end] for i in 1:n])) for j in 1:1000] .< r*sqrt(2/n)) 
#tests with alpha 0.01


n = 100
d= 2

#0.99 quantile vor mean of sum of squares of d-dimensional normals
r = sqrt(quantile(Distributions.Chisq(d), 0.99))


#  call 'sample(Wiener(d), t).X[:,i]' n times
function Wn(d, t, n, i)
 wn = zeros(d, n)
 for j in 1:n
   wn[:, j] = sample(Wiener(d), t ).X[:,i]
 end
 wn
end
  
function Bn(d, u, v, t, n, i)
 wn = zeros(d, n)
 for j in 1:n
   wn[:, j] = samplebridge(Wiener(d), u, v, t ).X[:,i]
 end
 wn
end
 
@Test.test size(sample(Wiener(d),  linspace(3.,4.,2)).X) == (d,2)
 
# the mean and variance of a Brownian motion at t=2 is 0 and 2
@Test.test norm(mean([ sample(Wiener(d),  linspace(0.,2.,5)).X[:,end] for i in 1:n])) < r*sqrt(2/n) #scale with std(W_2)/sqrt(n)
chiupper =  quantile(Distributions.Chisq(n), 0.995) #upper 0.005 percentile  
chilower = quantile(Distributions.Chisq(n), 0.005) #lower 0.005 percentile  
@Test.test chiupper >n*var( Wn(1, linspace(0.,2.,5), n, 5))/2 > chilower


# check that W(2) has the right covariance matrix
@Test.test (d-> norm(cov( Wn(d, linspace(0.,2.,5), n, 5)') - diagm(2ones(d)))*sqrt(n)/sqrt(d))(50) < 5. # did not figure that out exactly, should fail less then 99 %
@Test.test (d-> norm(cov( Wn(d,[0.,0.1, 0.3, 0.5, 2.], n, 5)') - diagm(2ones(d)))*sqrt(n)/sqrt(d))(50) < 5. # did not figure that out exactly, should fail less then 99 %


# a "deterministic" bridge with only start end endpoint
@Test.test (global t1 = norm(samplebridge(Wiener(2), [1.,2.], [5.,7.], linspace(1.,4.,2)).X - [1. 5.; 2. 7.])) < eps(10.)

# Check that the bridge has the right mean and Variance

# Covariance of a Brownian bridge from t_1 to t_2 (t_2-t)(s-t_1)/(t_2-t_1)
# here (3-0.5)*(0.5-0)/(3) = 0.4166666666666667


@Test.test norm(mean(Bn(2,[1.,2.], [5.,7.],[1.,1.1, 1.3, 1.5, 3.], n, 4),2) -  [1.,2.] - ([5.,7.] - [1.,2.])*.5/2) < sqrt(quantile(Distributions.Chisq(2), 0.99))*sqrt(0.416/n)

 # (3-0.5)*(0.5-0.1)/(2.9) = 0.3448275862068966
@Test.test quantile(Distributions.Chisq(n),  0.005) < n*var(Bn(1,[1.], [5.],[0.1, 0.3, 0.5, 3.], n, 3),2 )[1]/0.345  <  quantile(Distributions.Chisq(n), 0.995)

N = sample(PoissonProcess(1.), 0, 0., 10.1)
V = sample(Wiener(2), 0., linspace(0,1,5))
Vaug = augment( V, [0.3, 1., 10.])
dump(Vaug)
#julia> var(Bn(2,[1.,2.], [5.,7.],[0.1, 0.3, 0.5, 3.], 100000, 3),2 )
#2x1 Array{Float64,2}:
# 0.344571
# 0.343892



#W =  cumsum(dW,2) + (v .- sum(dW,2))*(t'.-t[1])/(t[end] - t[1]); VecPath(P, t, W)
#var([ sample(Wiener(2), zeros(2), linspace(0.,1.,3)).X[1,end] for i in 1:1000000]) \approx 1

