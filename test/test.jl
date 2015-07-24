import Anasol
using Base.Test

#test that the *_cf version matches with the *_c version when sourcestrength(t) = (inclosedinterval(t, t0, t1) ? 1. : 0.)

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
