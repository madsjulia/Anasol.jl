import Anasol
import LinearAlgebra
import DelimitedFiles
import Test

if !isdefined(Base, Symbol("@stderrcapture"))
	macro stderrcapture(block)
		quote
			if ccall(:jl_generating_output, Cint, ()) == 0
				errororiginal = Base.stderr;
				(errR, errW) = redirect_stderr();
				errorreader = @async read(errR, String);
				evalvalue = $(esc(block))
				redirect_stderr(errororiginal);
				close(errW);
				close(errR);
				return evalvalue
			end
		end
	end
end

@stderrcapture function contamination(wellx, welly, wellz, n, lambda, theta, vx, vy, vz, ax, ay, az, H, x, y, z, dx, dy, dz, f, t0, t1, t; anasolfunction=Anasol.long_bbb_ddd_iir_c)
	d = -theta * pi / 180
	xshift = wellx - x
	yshift = welly - y
	ztrans = wellz - z
	xtrans = xshift * cos(d) - yshift * sin(d)
	ytrans = xshift * sin(d) + yshift * cos(d)
	x01 = x02 = x03 = 0. # we transformed the coordinates so the source starts at the origin
	#sigma01 = sigma02 = sigma03 = 0.#point source
	sigma01 = dx
	sigma02 = dy
	sigma03 = dz
	v1 = vx
	v2 = vy
	v3 = vz
	speed = sqrt(vx * vx + vy * vy + vz * vz)
	sigma1 = sqrt(ax * speed * 2)
	sigma2 = sqrt(ay * speed * 2)
	sigma3 = sqrt(az * speed * 2)
	H1 = H2 = H3 = H
	xb1 = xb2 = xb3 = 0. # xb1 and xb2 will be ignored, xb3 should be set to 0 (reflecting boundary at z=0)
	anasolresult = anasolfunction([xtrans, ytrans, ztrans], t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1)
	return 1e6 * f * anasolresult / n
end

#a test using results that were verified against results from the C version of Mads/Anasol
@stderrcapture function testmadsc(anasolfunctionname)
	anasolfunction = eval(Meta.parse("Anasol.$anasolfunctionname"))
	resultsdir = string(dirname(Base.source_path()), "/goodresults")
	x, y, z = 1000, 1450, 0
	porosity = 0.1
	vx = 30.
	vz = vy = theta = lambda = 0.
	ax, ay, az = [70., 15., 0.3]
	H = 0.5
	dx, dy, dz = [250., 250., 1.]
	f = 50.
	t0, t1 = [5., 15.]
	wellx, welly, wellz = [823., 1499., 3.]
	ts = range(1.; stop=50., length=50)
	global results = Array{Float64}(undef, length(ts))
	for i = 1:length(ts)
		results[i] = contamination(wellx, welly, wellz, porosity, lambda, theta, vx, vy, vz, ax, ay, az, H, x, y, z, dx, dy, dz, f, t0, t1, ts[i]; anasolfunction=anasolfunction)
	end
	# DelimitedFiles.writedlm("$resultsdir/$anasolfunctionname.dat", results)
	goodresults = DelimitedFiles.readdlm("$resultsdir/$anasolfunctionname.dat")
	return LinearAlgebra.norm(results - goodresults)
end

@Test.testset "Anasol" begin
	x01, x02, x03 = 5., 3.14, 2.72
	global x0 = [x01, x02, x03]
	sigma01, sigma02, sigma03 = 1., 10., .1
	v1, v2, v3 = 2.5, 0.3, 0.01
	v = [v1, v2, v3]
	sigma1, sigma2, sigma3 = 100., 10., 1.
	H1, H2, H3 = 0.5, 0.5, 0.5
	xb1, xb2, xb3 = 0., 0., 0.
	lambda = 0.01
	t0, t1 = 0.5, sqrt(2)
	sourcestrength(t) = (Anasol.inclosedinterval(t, t0, t1) ? 1. : 0.)

	ts = range(0; stop=2, length=100)
	for t in ts
		for i = 1:1000
			x = x0 + v * t + 10 * randn(length(x0))
			@Test.test Anasol.long_bbb_ddd_iir_cf(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1, sourcestrength) == Anasol.long_bbb_ddd_iir_c(x, t, x01, sigma01, v1, sigma1, H1, xb1, x02, sigma02, v2, sigma2, H2, xb2, x03, sigma03, v3, sigma3, H3, xb3, lambda, t0, t1)
		end
	end

	t1s = collect(2015:5:2030)
	global results = Array{Float64}(undef, length(t1s))
	for i = 1:100
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
		flux = 50000
		t0 = 1985
		wellx = 1250.
		welly = 0.
		wellz0 = 3.
		wellz1 = 3.
		for t = range(2016; stop=2035, length=20)
			for j = 1:length(t1s)
				t1 = t1s[j]
				results[j] = .5 * (contamination(wellx, welly, wellz0, n, lambda, theta, vx, vy, vz, ax, ay, az, H, x, y, z, dx, dy, dz, flux, t0, t1, t) +
					contamination(wellx, welly, wellz1, n, lambda, theta, vx, vy, vz, ax, ay, az, H, x, y, z, dx, dy, dz, flux, t0, t1, t))
			end
			for j = 1:length(t1s) - 1
				@Test.test results[j] <= results[j + 1]
			end
		end
	end

	anasolfunctionnames = ["long_bbb_ddd_iir_c", "long_bbb_bbb_iir_c"]
	for anasolfunctionname in anasolfunctionnames
		@Test.test isapprox(testmadsc(anasolfunctionname), 0.; atol=1e-10)
	end
	include("newtest.jl")
end

:passed
