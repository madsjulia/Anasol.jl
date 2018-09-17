import Anasol
import Test
import DelimitedFiles

x = [0., 0.]
t = 1.
x0 = [0., 0.]
sigma0 = [1., 1.]
v = [1., 0.01]
sigma = [1., 1.]
H = [NaN, .75]
xb = [NaN, 0.]
t0 = 0.
t1 = 1.
lambda = 0.1
@stderrcapture function packargs(args...)
	returnargs = Array{Float64}(undef, length(args[1]) * length(args))
	for i = 1:length(args)
		for j = 1:length(args[1])
			returnargs[i + (j - 1) * length(args)] = args[i][j]
		end
	end
	return returnargs
end
disps = Val{(:linear,:fractional)}
sources = Val{(:dispersed,:box)}
boundaries = Val{(:infinite,:reflecting)}
distributions = [Anasol.standardnormal, Anasol.standardnormal]
Anasol.kernel(x, t, x0, sigma0, v, sigma, H, xb, disps, sources, boundaries)
@Test.test Anasol.kernel(x, t, x0, sigma0, v, sigma, H, xb, disps, sources, boundaries) == Anasol.long_bf_db_ir(x, t, packargs(x0, sigma0, v, sigma, H, xb)...)
@Test.test Anasol.kernel(x, t, x0, sigma0, v, sigma, H, xb, disps, sources, boundaries, distributions) == Anasol.long_bf_db_ir(x, t, packargs(x0, sigma0, v, sigma, H, xb)...)
@Test.test Anasol.kernel_c(x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, disps, sources, boundaries, distributions) == Anasol.long_bf_db_ir_c(x, t, packargs(x0, sigma0, v, sigma, H, xb)..., lambda, t0, t1)
@Test.test Anasol.kernel_cf(x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, t->exp(-t), disps, sources, boundaries, distributions) == Anasol.long_bf_db_ir_cf(x, t, packargs(x0, sigma0, v, sigma, H, xb)..., lambda, t0, t1, t->exp(-t))
N = 6
const dim = Val{2}
for i = 1:10 ^ N Anasol.innerkernel(dim, x, t, x0, sigma0, v, sigma, H, xb, disps, sources, boundaries, nothing) end
for i = 1:10 ^ N Anasol.innerkernel(dim, x, t, x0, sigma0, v, sigma, H, xb, disps, sources, boundaries, distributions) end
packedargs = packargs(x0, sigma0, v, sigma, H, xb)
for i = 1:10 ^ N Anasol.long_bf_db_ir(x, t, 0.0,1.0,1.0,1.0,NaN,NaN,0.0,1.0,0.01,1.0,0.75,0.0) end

for i = 1:10 ^ (N - 2) Anasol.kernel_c(x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, disps, sources, boundaries, distributions) end
for i = 1:10 ^ (N - 2) Anasol.long_bf_db_ir_c(x, t, 0.0,1.0,1.0,1.0,NaN,NaN,0.0,1.0,0.01,1.0,0.75,0.0, lambda, t0, t1) end

for i = 1:10 ^ (N - 2) Anasol.kernel_cf(x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, t->exp(-t), disps, sources, boundaries, distributions) end
for i = 1:10 ^ (N - 2) Anasol.long_bf_db_ir_cf(x, t, 0.0,1.0,1.0,1.0,NaN,NaN,0.0,1.0,0.01,1.0,0.75,0.0, lambda, t0, t1, t->exp(-t)) end