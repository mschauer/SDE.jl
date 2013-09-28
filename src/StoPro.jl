module StoPro
using Distributions
using NumericExtensions
#import getindex, setindex!
export MvWiener, PoissonProcess, DX, sample, samplebridge, law, at, augment

import Distributions.VariateForm
import Distributions.ValueSupport
typealias StateSpace ValueSupport


abstract IndexSet
type DiscreteTime <: IndexSet end
type ContinuousTime <: IndexSet end

# possible extensions:
#type GeneralStateSpace <: StateSpace end
#type GeneralIndexSet <: IndexSet end



# T time E type of elements of state vector
abstract StochasticProcess{T<:IndexSet, E<:StateSpace, V<:VariateForm}

# Systematic aliases
typealias DTChain StochasticProcess{DiscreteTime, Discrete, Univariate}
typealias CTChain StochasticProcess{ContinuousTime, Discrete, Univariate}
typealias DTProcess StochasticProcess{DiscreteTime, Continuous, Univariate}
typealias CTProcess StochasticProcess{ContinuousTime, Continuous, Univariate}
typealias CTVecProc StochasticProcess{ContinuousTime, Continuous, Multivariate}
typealias DTVecProc StochasticProcess{DiscreteTime, Continuous, Multivariate}

# Further aliases
typealias TimeSeries DTProcess
typealias MvTimeSeries DTVecProc



# generic sample path
abstract GenSamplePath




# type for path for which all jumps are a.s. 0 or 1
type CountingProcessPath  <:  GenSamplePath
	P :: CTChain
	t :: Array{Float64,1} #t contains the times of actual jumps of 1 = N[i+1]-N[i]  
				# and "observations" at times t[i] with 0 = N[i+1]-N[i]  
				# time intervals t[i] t[i+1] with N[i+1]-N[i] > 1 indicate unknown location of jumps of size 1, not as instantanious jump of size > 1
	N :: Array{Int64, 1}
end 

type VecPath  <:  GenSamplePath
	P :: CTVecProc
	t :: Array{Float64,1}
	X :: Array{Float64,2}

end 

getindex(V::CountingProcessPath, I) = (V.t[I], V.X[I])
getindex(V::VecPath, I) = (V.t[I], V.X[:,I])


type MvWiener <: CTVecProc
	d :: Int	
	 MvWiener(d::Integer) = d < 1 ? error("illegal dimension") : new(d)
end

type PoissonProcess <: CTChain
 	lambda :: Float64
	PoissonProcess(lambda :: Float64) = new(lambda)
end


function setindex!(V::VecPath,y, I) 
	V.t[I],V.X[:,I] = y 
#	getindex(V::VecPath, I) 
end




type SamplePathChain <:  GenSamplePath
	P :: CTChain
	t :: Array{Float64,1}
	X :: Array{Int64,1}
end

type DX <: GenSamplePath
	P :: CTVecProc 
	dt :: Array{Float64,1}
	dX :: Array{Float64,2}
end



dlin(t, n::Integer) = ones(n)*(t/n)





############ Poisson process ########################

function law(P::PoissonProcess, t)
 	Poisson(lambda*t)
end 

function at(C :: CountingProcessPath, t0::Float64)
	i = searchsortedlast(C.t, t0) #i == 0 throws bounds error 

	if C.t[i] != t0
		if (i == length(C.t)) || (C.N[i+1] - C.N[i] > 1)
	 		warning("possible jumps between $(C.t[i]) and $t0")
	 	end 
	end
		
	C.N[i]
end

function at(C :: CountingProcessPath, t::Vector) 
	[at(C, t0) for t0 in t]
end


function sample(P::PoissonProcess, n1, t1, t2)
	n2 = rand(Poisson(P.lambda*(t2-t1)))
	#there is no jump 
	samplebridge(P::PoissonProcess, t1, n1, t2, n2; jumpend = false)
end

function samplebridge(P::PoissonProcess, t1, n1, t2, n2; jumpend = true)
	
	dt = rand(Exponential(1./P.lambda), n2-n1 + 2)
	dt[1] = 0.0
	dt = dt/sum(dt)*(t2-t1)
	dt[1] = t1
	  
	
	N = collect(n1:n2+1)
	if (!jumpend)
	 	N[end] -= 1
	end
	CountingProcessPath(P, cumsum!(dt), N)
end

############ MvWiener process ########################


function augment(W :: VecPath, s) 
	P = W.P
	t = W.t #alias
	assert(issorted(s))
#	assert(min(t) <= min(s))  
	n = length(t)
	m = length(s)

	Wnew = VecPath(P, zeros(n + m), zeros(Float64, dim(P), n + m))
	
	j = 1 # traverse vector t 
	i = 1 #      and vector s
	while j <= n && i <= m #as long as there are entries in t and s
		assert(t[j] <= s[i]) #no sensible backward augmentation without starting distribution
		Wnew[j + i - 1] = W[j] #t[j] <= s[j], so append W[j] first

		if (i <= m && j == n) #something left in s and t exhausted
			Wnew.t[n + i : n + m] = s[i:m] #copy remaining time points

			#forward simulate from W.X[:, n] and append
			B = sample(P, W.X[:, n], Wnew.t[n + i-1 : n + m] ) 
			Wnew.X[:, n + i-1 : n + m] = B.X  

			break #s exhausted too

		elseif (i <= m && s[i] < t[j+1]   ) #entries left in s which should go before t[j+1]

			#which timepoints from s go inbetween t[j] and t[j+1]
			
			Wnew.t[j + i] = s[i]
			i2 = i
		  	while(i2 < m &&  s[i2+1] < t[j+1])  
		  		i2 += 1
	  			Wnew.t[j + i2] = s[i2]
		  	end
			Wnew.t[j + i2 + 1] = t[j+1]  
 
		  	# sample bridge from W[j] to W[j+1] using those timepoints
		  	B = samplebridge(P, W.X[:, j], W.X[:, j+1], Wnew.t[j + i - 1 : j + i2 + 1])
			Wnew.X[:,j + i - 1 : j + i2 + 1] = B.X
		
 			j += 1
		  	i = i2 + 1
		
		else 
			j += 1
		end
	  	
		
 	end
	Wnew
		 
end


function law(W::MvWiener, u, t)
	assert(dim(u) == d)
	IsoNormal(u, t)
end 

function dim(W::MvWiener)
	W.d
end
 
function sample(P::MvWiener, u, t)
	dt = [0.0, diff(t)]
	assert( min(dt) >= 0.0, "t linearly ordered")
	dW = randn(dim(P),length(dt)) .* sqrt(dt)'
	dW[:,1] = u
	VecPath(P, t, cumsum!(dW,2) )
end
function sample(P::MvWiener, t) 
	sample(P, zeros(dim(P)), t) 
end
 	

function samplebridge(P::MvWiener, u, v, t)
	dt = [0.0, diff(t)]
	assert( min(dt) >= 0.0, "t linearly ordered")
	dW = randn(dim(P),length(dt)) .* sqrt(dt)'
 	dW[:,1] = u 
 	dW = dW .+ (v .- sum(dW,2))*dt'/(t[end] - t[1])
  	VecPath(P, t, cumsum!(dW,2))
	
end

end
