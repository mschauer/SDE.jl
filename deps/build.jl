if OS_NAME == :Windows
	println("Windows: skipping binary deps, using julia fallback ...")
elseif OS_NAME == :Darwin
	println("Darwin: skipping binary deps, using julia fallback ...")
else
	println("Compiling sigma.c...")
	wd = pwd()
	cd(dirname(Base.source_path()))
	run(`gcc -fPIC -O3 -c sigma.c`)
	println("Linking libsigma...")
	run(`gcc -shared sigma.o -o libsigma.so`)
	cd(wd)
end
