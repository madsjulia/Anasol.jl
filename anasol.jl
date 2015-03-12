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
		kernelexpr = :(cdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau + .5 * $(symbol(string("sigma0", i)))) / $(symbol(string("sigmat", i)))) - cdf($(symbol(string("dist", i))), (x[$(i)] - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau - .5 * $(symbol(string("sigma0", i)))) / $(symbol(string("sigmat", i)))))
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

maxnumberofdimensions = 2
for n = 1:maxnumberofdimensions
	bigq = quote
		@nloops numberofdimensions j ii->1:length(boundarynames) begin
			@nloops numberofdimensions k ii->1:length(sourcenames) begin
				@nloops numberofdimensions i ii->1:length(dispersionnames) begin
					@nexprs numberofdimensions ii->(coreexprs_ii = coreexpression(dispersionnames[i_ii], dispersiontimedependence[i_ii], ii, sourcenames[k_ii], boundarynames[j_ii]))
					shortfunctionname = string((@ntuple numberofdimensions ii->dispersionnames[i_ii])..., "_", (@ntuple numberofdimensions ii->sourcenames[k_ii])..., "_", (@ntuple numberofdimensions ii->boundarynames[j_ii])...)
					q = quote
						$(symbol(string("long_", shortfunctionname)))(x::Vector,tau) = 1
					end
					returnexpr = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
					q.args[2].args[2].args[2] = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
					@nexprs numberofdimensions ii->q.args[2].args[2].args = [coreexprs_ii; q.args[2].args[2].args]
					for i = 1:numberofdimensions
						q.args[2].args[1].args = [q.args[2].args[1].args; symbol("x0$(i)"); symbol("sigma0$(i)"); symbol("v$(i)"); symbol("sigma$(i)"); symbol("H$(i)"); symbol("xb$(i)")]
					end
					eval(q)
					fullargs = copy(q.args[2].args[1].args)
					#now make a version that only includes the necessary arguments
					q.args[2].args[1].args = fullargs[1:3]
					q.args[2].args[1].args[1] = symbol(shortfunctionname)
					print(shortfunctionname, ": ")
					for i = 4:length(fullargs)
						argsymbol = fullargs[i]
						if Metatools.in(argsymbol, q.args[2].args[2])
							q.args[2].args[1].args = [q.args[2].args[1].args; argsymbol]
						end
					end
					eval(q)
					println(q.args[2].args[1].args)
				end
			end
		end
	end
	Metatools.replacesymbol!(bigq, :numberofdimensions, n)
	eval(bigq)
end

end
