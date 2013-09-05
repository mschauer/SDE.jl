function example1()
	println("euler")
	println(quvar(euler(0, 1, (t,x)-> -5x, (t,x) -> 1, diff(linspace(0., 1.,1000)))))

	dt = diff(linspace(0., 0.5, 100))
	dw = dW1(dt)

	println("Brownian motion with drift")
	println(" sint(4 .. dt) + ito (8 .. dw)")
	println((ito(4 .. dt) + ito(8 .. dw))[1:10])

	println(bracket(brown1(0, 2, 100))[1:10])

	#%discretization error of quadratic variation ?
	tic()
	println("Quadratic variation of B(0<t<2) = 2 \u2248 ", mean([quvar(brown1(0, 2, 5)) for i in 1:100000]))
	println("time:",toc())

	tic()
	dt2 = diff(linspace(0., 2,5))
	println("Quadratic variation of B(0<t<2) = 2 \u2248 ", mean([quvar(ito(dW1(dt2))) for i in 1:100000]))
	println("time:",toc())

	println("W_aug")
	adw = aug(dw,dt,10)
	println(adw[1:10])
	println("W(t) - W_aug(t):", sum(dw) - sum(adw))
	println("quvar(W, W_aug):", quvar(ito(1 .. dw)), " ", quvar(ito(1 .. adw)))

	dx = dW1(5,3)
	println(ito(dx .. dx)[end])
end
