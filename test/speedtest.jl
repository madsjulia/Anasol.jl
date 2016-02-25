import Anasol

function speedtest(N)
	x01, x02, x03 = 5., 3.14, 2.72
	x0 = [x01, x02, x03]
	sigma01, sigma02, sigma03 = 1., 10., .1
	sigma0 = [sigma01, sigma02, sigma03]
	v1, v2, v3 = 2.5, 0.3, 0.01
	v = [v1, v2, v3]
	sigma1, sigma2, sigma3 = 100., 10., 1.
	sigma = [sigma1, sigma2, sigma3]
	H1, H2, H3 = 0.5, 0.5, 0.5
	H = [H1, H2, H3]
	xb1, xb2, xb3 = 0., 0., 0.
	xb = [xb1, xb2, xb3]
	lambda = 0.01
	t0, t1 = 0.5, sqrt(2)
	ts = linspace(0, 2, 100)
	t = ts[1]
	x = x0 + v * t + 10 * randn(length(x0))
	xs = 5 * rand(length(x0), N)
	xs[1, :] += x01
	xs[2, :] += x02
	xs[3, :] += x03
	dispersions = Val{(:linear, :linear, :linear)}
	boxsources = Val{(:box, :box, :box)}
	gausssources = Val{(:dispersed, :dispersed, :dispersed)}
	boundaries = Val{(:infinite, :infinite, :infinite)}
	#@code_warntype Anasol.long_bbb_ddd_iir(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
	#@code_warntype Anasol.long_bbb_ddd_iir_ckernel(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
	fastruntime = @elapsed for t in ts
		for i = 1:N
			erf(xs[1, i])
			erf(xs[2, i])
			erf(xs[3, i])
		end
	end
	println("3 erf calls average run time: $(1000 * fastruntime / (length(ts) * N)) milliseconds")
	runtime = @elapsed for t in ts
		for i = 1:N
			Anasol.long_bbb_bbb_iii(xs[:, i], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
		end
	end
	println("Backwards compatible box source average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	println("factor: $(runtime / fastruntime)")
	runtime = @elapsed for t in ts
		for i = 1:N
			Anasol.innerkernel(Val{3}, xs[:, i], t, x0, sigma0, v, sigma, H, xb, dispersions, boxsources, boundaries, nothing)
		end
	end
	println("New box source average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	println("factor: $(runtime / fastruntime)")
	println()

	fastruntime = @elapsed for t in ts
		for i = 1:N
			exp(xs[1, i])
			exp(xs[2, i])
			exp(xs[3, i])
		end
	end
	println("3 exp calls average run time: $(1000 * fastruntime / (length(ts) * N)) milliseconds")
	runtime = @elapsed for t in ts
		for i = 1:N
			Anasol.long_bbb_ddd_iii(xs[:, i], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
		end
	end
	println("Backwards compatible Gaussian source average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	println("factor: $(runtime / fastruntime)")
	runtime = @elapsed for t in ts
		for i = 1:N
			Anasol.innerkernel(Val{3}, xs[:, i], t, x0, sigma0, v, sigma, H, xb, dispersions, gausssources, boundaries, nothing)
		end
	end
	println("New Gaussian source average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	println("factor: $(runtime / fastruntime)")
	println()
	fastruntime = runtime
	runtime = @elapsed for t in ts
		for i = 1:N
			Anasol.long_bbb_ddd_iii_c(xs[:, i], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1)
		end
	end
	println("Backwards compatible continuous release average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	println("factor: $(runtime / fastruntime) (compared to instantaneous source)")
	runtime = @elapsed for t in ts
		for i = 1:N
			Anasol.kernel_c(Val{3}, xs[:, i], t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, dispersions, gausssources, boundaries, nothing)
		end
	end
	println("New continuous release average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	println("factor: $(runtime / fastruntime)")
end

speedtest(round(Int, 1e4))
