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
abstract GenSamplePath{T<:IndexSet, E<:StateSpace, V<:VariateForm}



# Systematic aliases
typealias DTChain StochasticProcess{DiscreteTime, Discrete, Univariate}
typealias CTChain StochasticProcess{ContinuousTime, Discrete, Univariate}
typealias DTProcess StochasticProcess{DiscreteTime, Continuous, Univariate}
typealias CTProcess StochasticProcess{ContinuousTime, Continuous, Univariate}
typealias CTVecProc StochasticProcess{ContinuousTime, Continuous, Multivariate}
typealias DTVecProc StochasticProcess{DiscreteTime, Continuous, Multivariate}

typealias GenDTChainPath GenSamplePath{DiscreteTime, Discrete, Univariate}
typealias GenCTChainPath GenSamplePath{ContinuousTime, Discrete, Univariate}
typealias GenDTProcessPath GenSamplePath{DiscreteTime, Continuous, Univariate}
typealias GenCTProcessPath GenSamplePath{ContinuousTime, Continuous, Univariate}
typealias GenCTVecProcPath GenSamplePath{ContinuousTime, Continuous, Multivariate}
typealias GenDTVecProcPath GenSamplePath{DiscreteTime, Continuous, Multivariate}

# Further aliases
typealias TimeSeries DTProcess
typealias MvTimeSeries DTVecProc



# helper functions


#reversible diff
diff1(t) = [t[1], diff(t)]
#differential linspace
dlin(t, n::Integer) = ones(n)*(t/n)




# type for path for which all jumps are a.s. 0 or 1
type CountingProcessPath  <:  GenCTChainPath
	t :: Array{Float64,1} #t contains the times of actual jumps of 1 = N[i+1]-N[i]  
				# and "observations" at times t[i] with 0 = N[i+1]-N[i]  
				# time intervals t[i] t[i+1] with N[i+1]-N[i] > 1 indicate unknown location of jumps of size 1, not as instantanious jump of size > 1
	N :: Array{Int64, 1}
end 

type CTVecPath  <:  GenCTVecProcPath
	t :: Array{Float64,1}
	X :: Array{Float64,2}
end 

type SubCTVecPath  <:  GenCTVecProcPath
	t :: SubArray{Float64,1}
	X :: SubArray{Float64,2}
end 



getindex(V::GenCTChainPath, I) = (V.t[I], V.X[I])
getindex(V::GenCTVecProcPath, I) = (V.t[I], V.X[:,I])

function setindex!(V::GenCTVecPath,y, I) 
	V.t[I],V.X[:,I] = y 
end
function slice(V::GenCTVecPath, I) 
	SubCTVecPath(slice(V.t,I), slice(V.X, 1:dim(V.P), I)) 
end





### Processes 

type CorrWiener{Cov<:AbstractPDMat}  <: CTVecProc
	μ::Vector{Float64}
	Σ::Cov	
	d::Int	
	function MvWiener(μ::Vector{Float64}, Σ::Cov) 
		 d = length(μ)
		 dim(Σ) == d || throw(ArgumentError("The dimensions of μ and Σ are inconsistent."))
		 new(μ, Σ, d)
	end
end


type MvWiener <: CTVecProc
	d :: Int	
	MvWiener(d::Integer) = d < 1 ? error("illegal dimension") : new(d)
end


type MvWienerBridge <: CTVecProc
	d :: Int
	T :: Float64	
	v :: Vector{Float64}
	MvWienerBridge(T, v) = new(length(v), T, v)
end


type PoissonProcess <: CTChain
 	lambda :: Float64
	PoissonProcess(lambda :: Float64) = new(lambda)
end



type DX <: GenSamplePath
	P :: CTVecProc 
	dt :: Array{Float64,1}
	dX :: Array{Float64,2} #invariant: dX[1] = X[1]
	DX(P, dt, dX) = (assert( min(dt) >= 0.0, "negative stepsize");	new(P, dt, dX))
end

# sample white noise

function DW(u::Vector, t)
	d = length(u)
	dt = diff1(t)
	dW = randn(d,length(dt)) .* sqrt(dt)'
	dW[:,1]= u
	DX(MvWiener(d), dt, dW)
end

function DW!(X::Array{Float64,2}, u::Vector, t)
	dt = diff1(t)
	randn(X) .* sqrt(dt)'
	X[:,1]= u
	X
end





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
	  
	
	N = collect(int64(n1:n2+1))
	if (!jumpend)
	 	N[end] -= 1
	end
	CountingProcessPath(P, cumsum!(dt), N)
end

############ MvWiener process ########################



function law(W::MvWiener, u, t)
	assert(dim(u) == d)
	IsoNormal(u, t)
end 

function dim(W::MvWiener)
	W.d
end


function transform_noise!(P::MvWiener, V::SubVecPath)
	V.P = P
	cumsum!(V.dt)
	cumsum!(V.dX,2)
	V
end


function transform_noise!(P::MvWienerBridge, V :: SubVecPath)
	V.P = P
	
 	dW = dW .+ (P.v .- sum(dW,2))*dt'/P.T
	cumsum!(V.dt)
	assert(V.d
	cumsum!(V.dX,2)
	

	CTVecPath(P, cumsum(V.dt), cumsum!(V.dX,2) )
end


function sample!(P::MvWiener, u::Vector, V::VecPath)
	DW!(V.X, u, V.t)
	cumsum!(V.X,2)
	V	
end
 
function sample(P::MvWiener, u, t)
	dW = DW(u, t)
	CTVecPath(P, t, cumsum!(dW.dX,2) )
end

function sample(P::MvWiener, t) 
	sample(P, zeros(dim(P)), t) 
end
 	
function samplebridge(P::MvWiener, u, v, t)
	dW = DW(u, t)
	transform_noise!(MvWienerBridge(t[end]-t[1], v), dW)
  	CTVecPath(P, t, dW.dX)
end
function samplebridge!(P::MvWiener, u, v, V::SubVecPath)
	DW!(V.X, u, V.t)
	transform_noise!(MvWienerBridge(t[end]-t[1], v), dW)
  	CTVecPath(P, t, dW.dX)
end



function augment(P::GenCTProcess, W :: GenCTVecPath, s) 
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
			#B = sample(P, W.X[:, n], Wnew.t[n + i-1 : n + m] ) 
			#Wnew.X[:, n + i-1 : n + m] = B.X  
			sample!(P,  W.X[:, n], slice(Wnew, n + i-1 : n + m))
			
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
	#	  	B = samplebridge(P, W.X[:, j], W.X[:, j+1], Wnew.t[j + i - 1 : j + i2 + 1])
	#		Wnew.X[:,j + i - 1 : j + i2 + 1] = B.X
			samplebridge!(P,  W.X[:, n], W.X[:, j+1], slice(Wnew,j + i - 1 : j + i2 + 1))
		
 			j += 1
		  	i = i2 + 1
		
		else 
			j += 1
		end
	  	
		
 	end
	Wnew
		 
end


end
