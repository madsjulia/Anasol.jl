import Anasol

function speedtest(N)
	x01, x02, x03 = 5., 3.14, 2.72
	x0 = [x01, x02, x03]
	sigma01, sigma02, sigma03 = 1., 10., .1
	v1, v2, v3 = 2.5, 0.3, 0.01
	v = [v1, v2, v3]
	sigma1, sigma2, sigma3 = 100., 10., 1.
	H1, H2, H3 = 0.5, 0.5, 0.5
	xb1, xb2, xb3 = 0., 0., 0.
	lambda = 0.01
	t0, t1 = 0.5, sqrt(2)
	ts = linspace(0, 2, 100)
	t = ts[1]
	x = x0 + v * t + 10 * randn(length(x0))
	xs = 5 * rand(length(x0), N)
	xs[1, :] += x01
	xs[2, :] += x02
	xs[3, :] += x03
	#@code_warntype Anasol.long_bbb_ddd_iir(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
	#@code_warntype Anasol.long_bbb_ddd_iir_ckernel(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
	runtime = @elapsed @time for t in ts
		for i = 1:N
			erf(xs[1, i])
			erf(xs[2, i])
			erf(xs[3, i])
		end
	end
	println("3 erf calls average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	runtime = @elapsed @time for t in ts
		for i = 1:N
			Anasol.long_bbb_bbb_iii(xs[:, i], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
		end
	end
	println("Box source average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	runtime = @elapsed @time for t in ts
		for i = 1:N
			exp(xs[1, i])
			exp(xs[2, i])
			exp(xs[3, i])
		end
	end
	println("3 exp calls average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	runtime = @elapsed @time for t in ts
		for i = 1:N
			Anasol.long_bbb_ddd_iii(xs[:, i], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3)
		end
	end
	println("Gaussian source average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
	runtime = @elapsed @time for t in ts
		for i = 1:N
			Anasol.long_bbb_ddd_iii_c(xs[:, i], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1)
		end
	end
	println("Continuous release average run time: $(1000 * runtime / (length(ts) * N)) milliseconds")
end

speedtest(round(Int, 1e4))
