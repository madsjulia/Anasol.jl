__precompile__()

"""
MADS: Model Analysis & Decision Support in Julia (Mads.jl v1.0) 2016

http://mads.lanl.gov
http://madsjulia.lanl.gov
http://gitlab.com/mads/Mads.jl

Licensing: GPLv3: http://www.gnu.org/licenses/gpl-3.0.html

Copyright 2016.  Los Alamos National Security, LLC.  All rights reserved.

This material was produced under U.S. Government contract DE-AC52-06NA25396 for
Los Alamos National Laboratory, which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The Government is granted for itself and others acting on its
behalf a paid-up, nonexclusive, irrevocable worldwide license in this material to reproduce,
prepare derivative works, and perform publicly and display publicly. Beginning five (5) years after
--------------- November 17, 2015, ----------------------------------------------------------------
subject to additional five-year worldwide renewals, the Government is granted for itself and
others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
material to reproduce, prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.

NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY, LLC,
NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR
PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

LA-CC-15-080; Copyright Number Assigned: C16008
"""
module Anasol

using Distributions
using Base.Cartesian
using MetaProgTools

const standardnormal = Distributions.Normal(0, 1)

dispersionnames = ["b", "f"]#b is form brownian motion, f is for fractional brownian motion
dispersiontimedependence = [(factor, s)->:(sqrt($factor * $(symbol(string("sigma0", s))) ^ 2 + $(symbol(string("sigma", s))) ^ 2 * tau)), (factor, s)->:(sqrt($factor * $(symbol(string("sigma0", s))) ^ 2 + $(symbol(string("sigma", s))) ^ 2 * tau ^ (2 * $(symbol(string("H", s))))))]
#distributions = ["b"=>:standardnormal, "f"=>:standardnormal]
distributions = Dict(zip(["b"; "f"], [:standardnormal; :standardnormal]))
sourcenames = ["d", "b"]#d is for distributed (e.g., gaussian or alpha stable), b is for box
boundarynames = ["i", "r"]#d is for infinite (no boundary), r is for reflecting
functionnames = []

"Create core expressions"
function coreexpression(dispersionname, dispersiontimedependence, i, sourcename, boundaryname)
	q = quote
		global $(distributions[dispersionname])
		const $(symbol(string("dist", i))) = $(distributions[dispersionname])
		$(symbol(string("sigmat", i))) = $(dispersiontimedependence(sourcename == "b" ? 0 : 1, i))
		$(symbol(string("retval", i))) = fillintheblank
	end
	if sourcename == "d"
		kernelexpr = :(pdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau) / $(symbol(string("sigmat", i)))) / $(symbol(string("sigmat", i))))
	elseif sourcename == "b"
		kernelexpr = :((cdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau + .5 * $(symbol(string("sigma0", i)))) / $(symbol(string("sigmat", i)))) - cdf($(symbol(string("dist", i))), (replaceme - $(symbol(string("x0", i))) - $(symbol(string("v", i))) * tau - .5 * $(symbol(string("sigma0", i)))) / $(symbol(string("sigmat", i))))) / $(symbol(string("sigma0", i))))
	else
		error("Unknown source type: $sourcename")
	end
	if boundaryname == "i"
		MetaProgTools.replacesymbol!(kernelexpr, :replaceme, :(x[$(i)]))
	elseif boundaryname == "r"
		kernelexpr2 = deepcopy(kernelexpr)
		MetaProgTools.replacesymbol!(kernelexpr, :replaceme, :(x[$(i)]))
		MetaProgTools.replacesymbol!(kernelexpr2, :replaceme, :(2 * $(symbol(string("xb", i))) - x[$(i)]))
		kernelexpr = Expr(:call, :+, kernelexpr, kernelexpr2)
	else
		error("Unknown boundary condition: $boundaryname")
	end
	MetaProgTools.replacesymbol!(q, :fillintheblank, kernelexpr)
	return q
end

function inclosedinterval(x, a, b)
	return x >= a && x <= b
end

shortfunctionnames = []
maxnumberofdimensions = 3
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
					q.args[2].args[2].args[2] = parse(mapreduce(i->" * retval$i", *, "return retval1", 2:numberofdimensions))
					@nexprs numberofdimensions ii->q.args[2].args[2].args = [coreexprs_ii; q.args[2].args[2].args]
					for i = 1:numberofdimensions
						q.args[2].args[1].args = [q.args[2].args[1].args; symbol("x0$(i)"); symbol("sigma0$(i)"); symbol("v$(i)"); symbol("sigma$(i)"); symbol("H$(i)"); symbol("xb$(i)")]
					end
					#=
					if shortfunctionname == "bbb_bbb_iir"
						println("q:")
						println(q)
					end
					=#
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
					#=
					if shortfunctionname == "bbb_ddd_iir"
						println("qck:")
						println(qck)
					end
					=#
					eval(qck)#evaluate the kernel function definition
					#now make a function that integrates the kernel
					qc = quote
						function $(symbol(string("long_", shortfunctionname, "_c")))(thiswillbereplaced)#this function defines the continuous release function
							#these if statements let quadgk know where the discontinuities are
							if t - t0 <= 0
								#we are before the source started...return 0
								return 0.
							elseif t - t1 <= 0 && inclosedinterval(t - t0, 0, t)
								#the source has started, but not turned off...we don't need to integrate over (t-t0,t) because nothing has been in the system that long
								return quadgk(tau::Float64->$(symbol(string("long_", shortfunctionname, "_ckernel")))($([continuousreleaseargs[1:end]...; symbol("t")]...)), 0, t - t0; reltol=1e-7, abstol=1e-4)[1]
							elseif 0 <= t - t1 && t - t0 <= t
								#the source has started and turned off...we don't need to integrate over (0, t-t1) because everything has been in the system at least that long and we don't need to integrate over (t-t0,t) because nothing has been in the system that long
								return quadgk(tau::Float64->$(symbol(string("long_", shortfunctionname, "_ckernel")))($([continuousreleaseargs[1:end]...; symbol("t")]...)), t - t1, t - t0; reltol=1e-7, abstol=1e-4)[1]
							elseif inclosedinterval(t - t1, 0, t) && t - t0 >= t
								error("t0 is less than zero, but the code assumes that t0>=0")
							else
								error("outside of ifelses: [t, t0, t1] = [$t, $t0, $t1]")
							end
						end
					end
					continuousreleaseargs[2] = symbol("t")
					qc.args[2].args[1].args = [qc.args[2].args[1].args[1]; continuousreleaseargs[1:end]...]#give it the correct set of arguments
					#=
					if shortfunctionname == "bbb_ddd_iir"
						println("qc:")
						println(qc)
					end
					=#
					eval(qc)
					continuousreleaseargs[2] = symbol("tau")
					qcf = quote
						function $(symbol(string("long_", shortfunctionname, "_cf")))(thiswillbereplaced)#this function defines the continuous release function
							#these if statements let quadgk know where the discontinuities are
							if t - t0 <= 0
								#we are before the source started...return 0
								return 0.
							elseif t - t1 <= 0 && inclosedinterval(t - t0, 0, t)
								#the source has started, but not turned off...we don't need to integrate over (t-t0,t) because nothing has been in the system that long
								return quadgk(tau::Float64->sourcestrength(t - tau) * $(symbol(string("long_", shortfunctionname, "_ckernel")))($([continuousreleaseargs[1:end]...; symbol("t")]...)), 0, t - t0; reltol=1e-7, abstol=1e-4)[1]
							elseif 0 <= t - t1 && t - t0 <= t
								#the source has started and turned off...we don't need to integrate over (0, t-t1) because everything has been in the system at least that long and we don't need to integrate over (t-t0,t) because nothing has been in the system that long
								return quadgk(tau::Float64->sourcestrength(t - tau) * $(symbol(string("long_", shortfunctionname, "_ckernel")))($([continuousreleaseargs[1:end]...; symbol("t")]...)), t - t1, t - t0; reltol=1e-7, abstol=1e-4)[1]
							elseif inclosedinterval(t - t1, 0, t) && t - t0 >= t
								error("t0 is less than zero, but the code assumes that t0>=0")
							else
								error("outside of ifelses: [t, t0, t1] = [$t, $t0, $t1]")
							end
						end
					end
					continuousreleaseargs[2] = symbol("t")
					qcf.args[2].args[1].args = [qcf.args[2].args[1].args[1]; continuousreleaseargs[1:end]...; :(sourcestrength::Function)]#give it the correct set of arguments
					#now make a version that takes a function for the time-dependence of the source
					#=
					if shortfunctionname == "bbb_ddd_iir"
						println("qck:")
						println(qck)
						println("qcf:")
						println(qcf)
						println(qcf.args[2].args[1].args)
					end
					=#
					eval(qcf)
					#now make a version that only includes the necessary arguments
					fullargs = copy(q.args[2].args[1].args)
					q.args[2].args[1].args = fullargs[1:3]
					q.args[2].args[1].args[1] = symbol(shortfunctionname)
					for i = 4:length(fullargs)
						argsymbol = fullargs[i]
						if MetaProgTools.in(argsymbol, q.args[2].args[2])
							q.args[2].args[1].args = [q.args[2].args[1].args; argsymbol]
						end
					end
					#=
					if shortfunctionname == "bbb_ddd_iir"
						println("qgoodargs:")
						println(q)
					end
					=#
					eval(q)#make the function with only the necessary arguments
				end
			end
		end
	end
	MetaProgTools.replacesymbol!(bigq, :numberofdimensions, n)
	eval(bigq)
end

end
