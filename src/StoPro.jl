module StoPro
using Distributions
using NumericExtensions
#import getindex, setindex!
export MvWiener, PoissonProcess, diff1, sample!, sample, samplebridge, law, at, augment, DW, euler

import Distributions.VariateForm
import Distributions.ValueSupport
import Base.sub
import NumericExtensions.cumsum!
typealias StateSpace ValueSupport

# possible extensions:
#type GeneralStateSpace <: StateSpace end
#type GeneralIndexSet <: IndexSet end



# T time E type of elements of state vector
abstract StochasticProcess{E<:StateSpace, V<:VariateForm}
abstract GenSamplePath{E<:StateSpace, V<:VariateForm}



# Systematic aliases
typealias Chain StochasticProcess{ Discrete, Univariate}
typealias RealProcess StochasticProcess{ Continuous, Univariate}
typealias VecProc StochasticProcess{ Continuous, Multivariate}

typealias GenChainPath GenSamplePath{ Discrete, Univariate}
typealias GenRealProcessPath GenSamplePath{ Continuous, Univariate}
typealias GenVecProcPath GenSamplePath{ Continuous, Multivariate}

 


# helper functions


#reversible diff
diff1(a::StridedVector) = [i > 1 ? a[i] - a[i-1] : a[1] for i=1:length(a) ]
diff0(a::StridedVector) = [i > 1 ? a[i] - a[i-1] : zero(a[1]) for i=1:length(a) ]


#differential linspace
function dlin(t, T, n::Integer) 
	dt = ones(n)*(t/(n-1))
	dt[1] = t
	dt
end



# type for path for which all jumps are a.s. 0 or 1
type CountingProcessPath  <:  GenChainPath
	t :: Array{Float64,1} #t contains the times of actual jumps of 1 = N[i+1]-N[i]  
				# and "observations" at times t[i] with 0 = N[i+1]-N[i]  
				# time intervals t[i] t[i+1] with N[i+1]-N[i] > 1 indicate unknown location of jumps of size 1, not as instantanious jump of size > 1
	N :: Array{Int64, 1}
end 

type VecProcPath  <:  GenVecProcPath
	t :: Array{Float64,1}
	X :: Array{Float64,2}
end 

type SubVecProcPath  <:  GenVecProcPath
	t :: SubArray{Float64,1}
	X :: SubArray{Float64,2}
end 



getindex(V::GenChainPath, I) = (V.t[I], V.X[I])
getindex(V::GenVecProcPath, I) = (V.t[I], V.X[:,I])

function setindex!(V::GenVecProcPath,y, I) 
	V.t[I],V.X[:,I] = y 
end
function sub(V::GenVecProcPath, I) 
	SubVecProcPath(sub(V.t,I), sub(V.X, 1:size(V.X,1), I)) 
end




# crutch for missing inplace functions

function randn!{T,A}(X::SubArray{T,2,A})
	s = size(X)
	X[:,:] = randn(s...)
end

function cumsum!{T,A}(V::SubArray{T,1,A})
	V[:] = cumsum(V)
end

function cumsum!{T,A}(V::SubArray{T,2,A}, d)
	V[:,:] = cumsum(V,d)
end


### Processes 

type MvWhiteNoise <: VecProc
	d::Int	
	MvWhiteNoise(d::Integer) = d < 1 ? error("illegal dimension") : new(d)
end

type MvWiener <: VecProc
	d :: Int	
	MvWiener(d::Integer) = d < 1 ? error("illegal dimension") : new(d)
end


type CorrWiener{Cov<:AbstractPDMat}  <: VecProc
	μ::Vector{Float64}
	Σ::Cov	
	d::Int	
	function MvWiener(μ::Vector{Float64}, Σ::Cov) 
		 d = length(μ)
		 dim(Σ) == d || throw(ArgumentError("The dimensions of μ and Σ are inconsistent."))
		 new(μ, Σ, d)
	end
end

type MvWienerBridge <: VecProc
	d :: Int
	T :: Float64	
	v :: Vector{Float64}
	MvWienerBridge(T, v) = new(length(v), T, v)
end

type PoissonProcess <: Chain
 	lambda :: Float64
	PoissonProcess(lambda :: Float64) = new(lambda)
end




############ White Noise ########################

function DW(u::Vector, t)
	d = length(u)
	dt = diff0(t)
 	dW = randn(d,length(dt)) .* sqrt(dt)'
	dW[:,1]= u
	VecProcPath(t, dW)
end

function DW!(X::StridedArray{Float64,2}, u::Vector, t)
	dt = diff0(t)
 	randn!(X)
	X[:,:] = X .* sqrt(dt)'
	X[:,1]= u
	X
end

function DW!(V::GenVecProcPath,  u::Vector)
	dt = diff0(V.t)
 	randn(X) .* sqrt(dt)'
	V.X[:,1]= u
	V
end


function sample!(P::VecProc, u::Vector, V::GenVecProcPath)
	DW!(V.X, u, V.t)
	transformdw!(P, V)
end

 	
function sample(P::VecProc, u::Vector, t::Vector)
	dW = DW(u, t)
	transformdw!(P, dW)
end

function sample(P::VecProc, t) 
	sample(P, zeros(dim(P)), t) 
end
 	

function dim(V::VecProc)
	V.d
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
	CountingProcessPath(cumsum!(dt), N)
end

############ MvWiener process ########################



function law(W::MvWiener, u, t)
	assert(dim(u) == d)
	IsoNormal(u, t)
end 


function transformdw!(P::MvWiener, V::GenVecProcPath)
	cumsum!(V.X,2)
	V
end


function transformdw!(P::MvWienerBridge, V :: GenVecProcPath) #here	u = V.X[:, 1]
	assert(P.T == V.t[end])
 	V.X[:] = V.X .+ (P.v - sum(V.X,2))*diff0(V.t)'/(P.T - V.t[1])
	cumsum!(V.X,2)
	V
end

 

function samplebridge(P::MvWiener, u, v, t)
	dW = DW(u, t)
	transformdw!(MvWienerBridge(t[end], v), dW)
  	dW
end
function samplebridge!(P::MvWiener, u, v, V::GenVecProcPath)
	DW!(V.X, u, V.t)
	transformdw!(MvWienerBridge(V.t[end], v),V)
	V
end

############ Retrospective subsampling of time  ########################

function augment(P::VecProc, W :: GenVecProcPath, s) 
	t = W.t 
	assert(issorted(s))
	n = length(t)
	m = length(s)

	Wnew = VecProcPath(zeros(n + m), zeros(Float64, dim(P), n + m))
	
	j = 1 # traverse vector t 
	i = 1 #      and vector s
	while j <= n && i <= m #as long as there are entries in t and s
		assert(t[j] <= s[i]) #no sensible backward augmentation without starting distribution
		Wnew[j + i - 1] = W[j] #t[j] <= s[j], so append W[j] first

		if (i <= m && j == n) #something left in s and t exhausted
			Wnew.t[n + i : n + m] = s[i:m] #copy remaining time points

			#forward simulate from W.X[:, n] and append
			sample!(P,  W.X[:, n], sub(Wnew, n + i-1 : n + m))
			
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
			samplebridge!(P,  W.X[:, n], W.X[:, j+1], sub(Wnew,j + i - 1 : j + i2 + 1))
		
 			j += 1
		  	i = i2 + 1
		
		else 
			j += 1
		end
	  	
		
 	end
	Wnew
		 
end



#%  .. function:: euler(u, b(t,i,x), sigma(t,i,x), V :: GenVecProcPath)
#%  
#%  	Multivariate euler scheme, starting in u.
#% 	``dw`` -- Wiener differential with ``n`` values in the 
#%  	the interval ``[t0,sum(dt)]`` sampled at timepoints ``t0+Dt[1], t0 + Dt[1] + Dt[2], ...``
#%	``b, sigma`` -- drift and diffusion coefficient. 
#%	To allow for
#%  	
 		
function euler(u, b, sigma, V :: GenVecProcPath)
	t = V.t 
	dW = V.X
	
	N = (size(dW))[end]
	
	X = zeros(length(u), N)

	y = copy(u)
	for i in 1:N-1
		X[:,i] = y
 		y[:] = y .+  b(t,i, y)*(t[i+1]-t[i]) .+ sigma(t,i, y)*dW[:,i+1]
	
	end
	X[:,N] = y
	
	VecProcPath(t, X)
end




end

