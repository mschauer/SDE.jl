k = 0.5
Z1(X) = sum(X[1:end-1].^2.*diff(atan(X)))
Z2(X) = sum(X[1:end-1].^2.*diff(sin(X)))
Z3(X) = sum(X[1:end-1].^2.*0.01)

Z1() = Z1(cumsum(0.1*randn(100)))
Z2() = Z2(cumsum(0.1*randn(100)))
Z3() = Z3(cumsum(0.1*randn(100)))

var([(k*Z1()) for i in 1:10000])
var([(k*Z2()) for i in 1:10000])
var([(k*Z3()) for i in 1:10000])

mean([exp(k*Z1()) for i in 1:5000])
mean([exp(k*Z1()) for i in 1:50000])
mean([exp(k*Z3()) for i in 1:5000])
mean([exp(k*Z3()) for i in 1:50000])
srand(3);mean([exp(k*Z2()) for i in 1:5000])
srand(3);mean([exp(k*Z2()) for i in 1:50000])
#srand(3);mean([exp(k*Z2()) for i in 1:500000])

var([exp(k*Z1()) for i in 1:5000])
var([exp(k*Z1()) for i in 1:50000])
var([exp(k*Z2()) for i in 1:5000])
var([exp(k*Z2()) for i in 1:50000])
var([exp(k*Z3()) for i in 1:5000])
var([exp(k*Z3()) for i in 1:50000])
