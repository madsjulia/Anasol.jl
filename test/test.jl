import Anasol
import Mads
using Base.Test

#test that the *_cf version matches with the *_c version when sourcestrength(t) = (inclosedinterval(t, t0, t1) ? 1. : 0.)
function testcf()
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
	sourcestrength(t) = (Anasol.inclosedinterval(t, t0, t1) ? 1. : 0.)

	ts = linspace(0, 2, 100)
	for t in ts
		for i = 1:1000
			x = x0 + v * t + 10 * randn(length(x0))
			@test Anasol.long_bbb_ddd_iir_cf(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1, sourcestrength) == Anasol.long_bbb_ddd_iir_c(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1)
		end
	end
end

#tests that if we run the source longer, the concentrations increase
function testmonotone(N)
	t1s = collect(2015:5:2030)
	results = Array(Float64, length(t1s))
	for i = 1:N
		n = 0.1
		lambda = 0.
		theta = 0.
		vx = 30. + 10 * rand() - 5
		vy = 0.
		vz = 0.
		ax = .5 * 175. + 105 * rand() - 52.5
		ay = 15.
		az = 0.3
		H = 0.5
		x = 0.
		y = 100.
		z = 0.
		dx = 250.
		dy = 100.
		dz = 1.
		f = 50000
		t0 = 1985
		wellx = 1250.
		welly = 0.
		wellz0 = 3.
		wellz1 = 3.
		usefff = false
		for t = linspace(2016, 2035, 20)
			for j = 1:length(t1s)
				t1 = t1s[j]
				results[j] = .5 * (Mads.contamination(wellx, welly, wellz0, n, lambda, theta, vx, vy, vz, ax, ay, az, H, x, y, z, dx, dy, dz, f, t0, t1, t, usefff) + Mads.contamination(wellx, welly, wellz1, n, lambda, theta, vx, vy, vz, ax, ay, az, H, x, y, z, dx, dy, dz, f, t0, t1, t, usefff))
			end
			for j = 1:length(t1s) - 1
				@test results[j] <= results[j + 1]
			end
		end
	end
end

testcf()
testmonotone(100)
