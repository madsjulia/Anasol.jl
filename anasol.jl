module Anasol

using Distributions
using Base.Cartesian
using Metatools

const standardnormal = Distributions.Normal(0, 1)

dispersionnames = ["b", "f"]#b is form brownian motion, f is for fractional brownian motion
dispersiontimedependence = [s->:(sqrt($(symbol(string("sigma0", s))) + $(symbol(string("sigma", s))) * tau)), s->:(sqrt($(symbol(string("sigma0", s))) + $(symbol(string("sigma", s))) * tau ^ (2 * $(symbol(string("H", s))))))]
distributions = ["b"=>:standardnormal, "f"=>:standardnormal]
sourcenames = ["d", "b"]#d is for distributed (e.g., gaussian or alpha stable), b is for box
boundarynames = ["i", "r"]#d is for infinite (no boundary), r is for reflecting
functionnames = []

function coreexpression(dispersionname, dispersiontimedependence, i, sourcename, boundaryname)
	q = quote
		global $(distributions[dispersionname])
		const $(symbol(string("dist", i))) = $(distributions[dispersionname])
		$(symbol(string("sigmat", i))) = $(dispersiontimedependence(i))
		$(symbol(string("retval", i))) = fillintheblank
	end
	if sourcename == "d"
		kernelexpr = :(pdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau) / $(symbol(string("sigmat", i)))) / $(symbol(string("sigmat", i))))
	elseif sourcename == "b"
		kernelexpr = :(cdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau + .5 * $(symbol(string("sigma0", i)))) / $(symbol(string("sigmat", i)))) - cdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau - .5 * $(symbol(string("sigma0", i)))) / $(symbol(string("sigmat", i)))))
	else
		error("Unknown source type: $sourcename")
	end
	if boundaryname == "i"
		Metatools.replacesymbol!(kernelexpr, :replaceme, :(x[$(i)]))
	elseif boundaryname == "r"
		kernelexpr2 = deepcopy(kernelexpr)
		Metatools.replacesymbol!(kernelexpr, :replaceme, :(x[$(i)]))
		Metatools.replacesymbol!(kernelexpr2, :replaceme, :($(symbol(string("xb", i))) - x[$(i)]))
		kernelexpr = Expr(:call, :+, kernelexpr, kernelexpr2)
	else
		error("Unknown boundary condition: $boundaryname")
	end
	Metatools.replacesymbol!(q, :fillintheblank, kernelexpr)
	return q
end

function inclosedinterval(x, a, b)
	return x >= a && x <= b
end

shortfunctionnames = []
maxnumberofdimensions = 2
for n = 1:maxnumberofdimensions
	bigq = quote
		@nloops numberofdimensions j ii->1:length(boundarynames) begin
			@nloops numberofdimensions k ii->1:length(sourcenames) begin
				@nloops numberofdimensions i ii->1:length(dispersionnames) begin
					@nexprs numberofdimensions ii->(coreexprs_ii = coreexpression(dispersionnames[i_ii], dispersiontimedependence[i_ii], ii, sourcenames[k_ii], boundarynames[j_ii]))
					shortfunctionname = string((@ntuple numberofdimensions ii->dispersionnames[i_ii])..., "_", (@ntuple numberofdimensions ii->sourcenames[k_ii])..., "_", (@ntuple numberofdimensions ii->boundarynames[j_ii])...)
					shortfunctionnames = [shortfunctionnames; shortfunctionname]
					q = quote
						$(symbol(string("long_", shortfunctionname)))(x::Vector,tau) = 1
					end
					returnexpr = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
					q.args[2].args[2].args[2] = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
					@nexprs numberofdimensions ii->q.args[2].args[2].args = [coreexprs_ii; q.args[2].args[2].args]
					for i = 1:numberofdimensions
						q.args[2].args[1].args = [q.args[2].args[1].args; symbol("x0$(i)"); symbol("sigma0$(i)"); symbol("v$(i)"); symbol("sigma$(i)"); symbol("H$(i)"); symbol("xb$(i)")]
					end
					println("q:")
					println(q)
					eval(q)#make the function with all possible arguments
					#now make a version that includes a continuously released source from t0 to t1
					continuousreleaseargs = [q.args[2].args[1].args[2:end]; symbol("lambda"); symbol("t0"); symbol("t1")]
					#start by making the kernel of the time integral
					qck = quote
						function $(symbol(string("long_", shortfunctionname, "_ckernel")))(thiswillbereplaced)#this function defines the kernel that the continuous release function integrates against
							if inclosedinterval(t - tau, t0, t1)
								retval = exp(-lambda * tau) * $(symbol(string("long_", shortfunctionname)))($((q.args[2].args[1].args[2:end])...))
								return retval
							else
								return 0.
							end
						end
					end
					qck.args[2].args[1].args = [qck.args[2].args[1].args[1]; continuousreleaseargs[1:end]...; symbol("t")]#give it the correct set of arguments
					println("qck:")
					println(qck)
					eval(qck)#evaluate the kernel function definition
					#now make a function that integrates the kernel
					qc = quote
						function $(symbol(string("long_", shortfunctionname, "_c")))(thiswillbereplaced)#this function defines the kernel that the continuous release function integrates against
							return quadgk(tau->$(symbol(string("long_", shortfunctionname, "_ckernel")))($([continuousreleaseargs[1:end]...; symbol("t")]...)), 0, t)[1]
						end
					end
					continuousreleaseargs[2] = symbol("t")
					qc.args[2].args[1].args = [qc.args[2].args[1].args[1]; continuousreleaseargs[1:end]...]#give it the correct set of arguments
					println("qc:")
					println(qc)
					eval(qc)
					#now make a version that only includes the necessary arguments
					fullargs = copy(q.args[2].args[1].args)
					q.args[2].args[1].args = fullargs[1:3]
					q.args[2].args[1].args[1] = symbol(shortfunctionname)
					for i = 4:length(fullargs)
						argsymbol = fullargs[i]
						if Metatools.in(argsymbol, q.args[2].args[2])
							q.args[2].args[1].args = [q.args[2].args[1].args; argsymbol]
						end
					end
					eval(q)#make the function with only the necessary arguments
				end
			end
		end
	end
	Metatools.replacesymbol!(bigq, :numberofdimensions, n)
	eval(bigq)
end

end
