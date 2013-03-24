function example1()
	println("euler")
	println(quvar(euler(0, 1, (t,x)-> -5x, (t,x) -> 1, dT(1,1000))))

	dt = dT(0.5, 100)
	dw = dW(dt)

	println("Brownian motion with drift")
	println(" sint(4 .. dt) + ito (8 .. dw)")
	println(sint(4 .. dt) + ito (8 .. dw))

	println(bracket(brown(0, 2, 100)))

	#%test: discretization error of quadratic variation
	tic()
	println("Quadratic variation of B(0<t<2) = 2 \u2248 ", mean([quvar(brown(0, 2, 5)) for i in 1:100000]))
	println("time:",toc())

	tic()
	dt2 = dT(2, 5)
	println("Quadratic variation of B(0<t<2) = 2 \u2248 ", mean([quvar(ito(dW(dt2))) for i in 1:100000]))
	println("time:",toc())

	println("W_aug")
	adw = aug(dw,dt,10)
	println(adw)
	println("W(t) - W_aug(t):", sum(dw) - sum(adw))
	println("quvar(W, W_aug):", quvar(ito(1 .. dw)), " ", quvar(ito(1 .. adw)))

	dx = dW(5,10,3)
	println(ito(dx .. dx)[end])
end
