module Anasol

using Distributions
using Base.Cartesian
using Metatools

const standardnormal = Distributions.Normal(0, 1)

dispersionnames = ["b", "f"]#b is form brownian motion, f is for fractional brownian motion
dispersiontimedependence = [s->:(sqrt(tau)), s->:(tau ^ $(symbol(string("H", s))))]
distributions = ["b"=>:standardnormal, "f"=>:standardnormal]
functionnames = []

function coreexpression(dispersionname, dispersiontimedependence, i)
	return quote
		global $(distributions[dispersionname])
		const $(symbol(string("dist", i))) = $(distributions[dispersionname])
		$(symbol(string("sigmat", i))) = $(symbol(string("sigma", i))) * $(dispersiontimedependence(i))
		$(symbol(string("retval", i))) = pdf($(symbol(string("dist", i))), (x[$(i)] - $(symbol(string("v", i))) * tau) / $(symbol(string("sigmat", i)))) / $(symbol(string("sigmat", i)))
	end
end

maxnumberofdimensions = 3
for n = 1:maxnumberofdimensions
	bigq = quote
		@nloops numberofdimensions i k->1:length(dispersionnames) begin
			@nexprs numberofdimensions j->(coreexprs_j = coreexpression(dispersionnames[i_j], dispersiontimedependence[i_j], j))
			functionname = string((@ntuple numberofdimensions k->dispersionnames[i_k])...)
			q = quote
				$(symbol(functionname))(x::Vector,tau) = 1
			end
			returnexpr = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
			q.args[2].args[2].args[2] = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
			@nexprs numberofdimensions j->q.args[2].args[2].args = [coreexprs_j; q.args[2].args[2].args]
			for ii = 1:numberofdimensions
				q.args[2].args[1].args = [q.args[2].args[1].args; symbol("v$(ii)"); symbol("sigma$(ii)"); symbol("H$(ii)")]
			end
			eval(q)
		end
	end
	Metatools.replacesymbolwithvalue!(bigq, :numberofdimensions, n)
	eval(bigq)
end

end
