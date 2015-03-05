module Anasol

using Distributions

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

for i = 1:length(dispersionnames)
	functionname = string(dispersionnames[i])
	functionnames = [functionnames; functionname]
	eval(quote
		function $(symbol(functionname))(x, tau, v1, sigma1, H1)
			$(coreexpression(dispersionnames[i], dispersiontimedependence[i], 1))
			return retval1
		end
	end)
end

for i = 1:length(dispersionnames)
	for j = 1:length(dispersionnames)
		functionname = string(dispersionnames[i], dispersionnames[j])
		functionnames = [functionnames; functionname]
		eval(quote
			function $(symbol(functionname))(x::Vector, tau, v1, sigma1, H1, v2, sigma2, H2)
				$(coreexpression(dispersionnames[i], dispersiontimedependence[i], 1))
				$(coreexpression(dispersionnames[j], dispersiontimedependence[j], 2))
				return retval1 * retval2
			end
		end)
	end
end

for i = 1:length(dispersionnames)
	for j = 1:length(dispersionnames)
		for k = 1:length(dispersionnames)
			functionname = string(dispersionnames[i], dispersionnames[j], dispersionnames[k])
			functionnames = [functionnames; functionname]
			q = quote
				function $(symbol(functionname))(x::Vector, tau, v1, sigma1, H1, v2, sigma2, H2, v3, sigma3, H3)
					$(coreexpression(dispersionnames[i], dispersiontimedependence[i], 1))
					$(coreexpression(dispersionnames[j], dispersiontimedependence[j], 2))
					$(coreexpression(dispersionnames[j], dispersiontimedependence[j], 3))
					return retval1 * retval2 * retval3
				end
			end
			eval(q)
		end
	end
end



end
