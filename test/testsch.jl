#include(Pkg.dir("SDE","src", "Schauder.jl"))
#using Schauder
using SDE.Schauder
using Base.Test

th1 = sin(linspace(0,10,2^5-1))
th2 = copy(th1)
pickup!(th2)
drop!(th2)
@test norm(th1 - th2) < 2.0*eps()

@test [level(ones(i)) for i in 1:15] ==  [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4]

@test all([vectoroflevels(5,0),  5, vectoroflevels(5,0)] .== vectoroflevels(6,0))

L = 4
@test [number(L,0), 2^L-1+(1:5)] == number(L, 5)

@test unique(diff(vectoroflevels(L,0)[sortperm(number(L,0))])) == [0,1]


@test_approx_eq norm(fe_transf(sin, 0., 2pi, 3).^2 + fe_transf(sin, pi/2, 2.5pi, 3).^2 -1.0) 0

@test norm(fe_transfB1(sin, 0, 2pi, L) - [fe_transf(sin, 0, 2pi, L), 0]) < 4eps()

@test norm(fe_transfB2(sin, 0, 2pi, L) - [fe_transf(sin, 0, 2pi, L),0,0]) < 4eps()

@test dropB1([0,2,0,1]) == [1,3,1,1]

@test dropB2([0.,2.,1.,1.,1.]) == [1., 1., 3., 2., 1.]

