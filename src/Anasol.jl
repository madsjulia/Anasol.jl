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

import MetaProgTools
import DocumentFunction
import Distributions
import QuadGK
import DelimitedFiles
import Compat
import Compat.string
using Base.Cartesian

include("newanasol.jl")

const standardnormal = Distributions.Normal(0, 1)

documentation = false

docarguments = """Arguments

- `t`: time to compute concentration
- `x`: spatial coordinates of the point to compute the concentration
- `x01`/`x02`/`x03`: contaminant source coordinates
- `sigma01`/`sigma02`/`sigma01`: contaminant source sizes (if a constrained source) or standard deviations (if a distributed source)
- `sourcestrength`: user-provided function defining time-dependent source strength
- `t0`/`t1`: contaminant release times (source is released  between `t0` and `t1`)
- `v1`/`v2`/`v3`: groundwater flow velocity components
- `sigma1`/`sigma2`/`sigma3`: groundwater flow dispersion components
- `lambda`: half-life contaminant decay
- `H1`/`H2`/`H3`: Hurst coefficients in the case of fractional Brownian dispersion
- `xb1`/`xb2`/`xb3`: locations of the domain boundaries

Returns:

- contaminant concentration at location `x` at time `t`
"""

dictarguments = Dict(
"t"=>"time to compute concentration",
"x"=>"spatial coordinates of the point to compute the concentration",
"x01"=>"`x` contaminant source coordinate",
"x02"=>"`y` contaminant source coordinate",
"x03"=>"`z` contaminant source coordinate",
"sigma01"=>"`x` contaminant source size (if a constrained source) or standard deviation (if a distributed source)",
"sigma02"=>"`y` contaminant source size (if a constrained source) or standard deviation (if a distributed source)",
"sigma03"=>"`z` contaminant source size (if a constrained source) or standard deviation (if a distributed source)",
"t0"=>"contaminant release start time",
"t1"=>"contaminant release end time",
"tau"=>"integration constant",
"v1"=>"`x` groundwater flow velocity",
"v2"=>"`y` groundwater flow velocity",
"v3"=>"`z` groundwater flow velocity",
"sigma1"=>"`x` groundwater flow dispersion",
"sigma2"=>"`y` groundwater flow dispersion",
"sigma3"=>"`z` groundwater flow dispersion",
"sourcestrength"=>"user-provided function defining time-dependent source strength",
"lambda"=>"half-life contaminant decay",
"H1"=>"`x` Hurst coefficient in the case of fractional Brownian dispersion",
"H2"=>"`y` Hurst coefficient in the case of fractional Brownian dispersion",
"H3"=>"`z` Hurst coefficient in the case of fractional Brownian dispersion",
"xb1"=>"`x` location of the domain boundary",
"xb2"=>"`y` location of the domain boundary",
"xb3"=>"`z` location of the domain boundary",
)

function inclosedinterval(x, a, b)
	return x >= a && x <= b
end

axisnames = ["X", "Y", "Z"]
dispersionnames = ["b", "f"] # b is form brownian motion, f is for fractional brownian motion
sourcenames = ["d", "b"] # d is for distributed (e.g., Gaussian or Levy alpha stable), b is for box
boundarynames = ["i", "r", "a"] # d is for infinite (no boundary), r is for reflecting

function getdispersions(dispersionnames)
	f(x) = x == "b" ? :linear : :fractional
	return :(Val{$(map(f, dispersionnames))})
end
function getsources(sourcenames)
	f(x) = x == "b" ? :box : :dispersed
	return :(Val{$(map(f, sourcenames))})
end
function getboundaries(boundarynames)
	f(x) = x == "i" ? :infinite : x == "r" ? :reflecting : :absorbing
	return :(Val{$(map(f, boundarynames))})
end

function getlongdispersions(dispersionnames)
	s = "Dispersion:"
	for i = 1:length(dispersionnames)
		s *= " " * axisnames[i] * " - "
		s *= dispersionnames[i] == "b" ? "brownian" : "fractional"
		s *= " ($(dispersionnames[i]))"
	end
	return s
end
function getlongsources(sourcenames)
	s = "Source:"
	for i = 1:length(sourcenames)
		s *= " " * axisnames[i] * " - "
		s *= sourcenames[i] == "b" ? "constrained" : "dispersed"
		s *= " ($(sourcenames[i]))"
	end
	return s
end
function getlongboundaries(boundarynames)
	s = "Boundary:"
	for i = 1:length(boundarynames)
		s *= " " * axisnames[i] * " - "
		s *= boundarynames[i] == "b" ? "infinite" : boundarynames[i] == "r" ? "reflecting" : "absorbing"
		s *= " ($(boundarynames[i]))"
	end
	return s
end

#the functions defined in this monstrosity of loops are for backwards compatibility
#use the functions defined in "newanasol.jl" instead
maxnumberofdimensions = 3
for n = 1:maxnumberofdimensions
	bigq = quote
		@nloops numberofdimensions j ii->1:length(boundarynames) begin
			@nloops numberofdimensions k ii->1:length(sourcenames) begin
				@nloops numberofdimensions i ii->1:length(dispersionnames) begin
					shortfunctionname = string((@ntuple numberofdimensions ii->dispersionnames[i_ii])..., "_", (@ntuple numberofdimensions ii->sourcenames[k_ii])..., "_", (@ntuple numberofdimensions ii->boundarynames[j_ii])...)
					q = quote
						$(Symbol(string("long_", shortfunctionname)))(x::Vector,tau) = 1
					end
					x0s = Meta.parse(string("[", join(map(i->"x0$i", 1:numberofdimensions), ",")..., "]"))
					sigma0s = Meta.parse(string("[", join(map(i->"sigma0$i", 1:numberofdimensions), ",")..., "]"))
					vs = Meta.parse(string("[", join(map(i->"v$i", 1:numberofdimensions), ",")..., "]"))
					sigmas = Meta.parse(string("[", join(map(i->"sigma$i", 1:numberofdimensions), ",")..., "]"))
					Hs = Meta.parse(string("[", join(map(i->"H$i", 1:numberofdimensions), ",")..., "]"))
					xbs = Meta.parse(string("[", join(map(i->"xb$i", 1:numberofdimensions), ",")..., "]"))
					dispersions = getdispersions(@ntuple numberofdimensions ii->dispersionnames[i_ii])
					docdispersions = getlongdispersions(@ntuple numberofdimensions ii->dispersionnames[i_ii])
					sources = getsources(@ntuple numberofdimensions ii->sourcenames[k_ii])
					docsources = getlongsources(@ntuple numberofdimensions ii->sourcenames[k_ii])
					boundaries = getboundaries(@ntuple numberofdimensions ii->boundarynames[j_ii])
					docboundaries = getlongboundaries(@ntuple numberofdimensions ii->boundarynames[j_ii])
					q.args[2].args[2].args[2] = :(innerkernel(Val{$numberofdimensions}, x, tau, $x0s, $sigma0s, $vs, $sigmas, $Hs, $xbs, $dispersions, $sources, $boundaries, nothing))
					for i = 1:numberofdimensions
						q.args[2].args[1].args = [q.args[2].args[1].args; Symbol("x0$(i)"); Symbol("sigma0$(i)"); Symbol("v$(i)"); Symbol("sigma$(i)"); Symbol("H$(i)"); Symbol("xb$(i)")]
					end
					eval(q)# make the function with all possible arguments
					# now make a version that includes a continuously released source from t0 to t1
					continuousreleaseargs = [q.args[2].args[1].args[2:end]; Symbol("lambda"); Symbol("t0"); Symbol("t1")]
					# start by making the kernel of the time integral
					qck = quote
						function $(Symbol(string("long_", shortfunctionname, "_ckernel")))(thiswillbereplaced) # this function defines the kernel that the continuous release function integrates against
							return cinnerkernel(Val{$numberofdimensions}, x, tau, $x0s, $sigma0s, $vs, $sigmas, $Hs, $xbs, lambda, t0, t1, t, $dispersions, $sources, $boundaries, nothing)
						end
					end
					qck.args[2].args[1].args = [qck.args[2].args[1].args[1]; continuousreleaseargs[1:end]...; Symbol("t")] # give it the correct set of arguments
					eval(qck) # evaluate the kernel function definition
					# now make a function that integrates the kernel
					qc = quote
						function $(Symbol(string("long_", shortfunctionname, "_c")))(thiswillbereplaced) # this function defines the continuous release function
							return kernel_c(x, t, $x0s, $sigma0s, $vs, $sigmas, $Hs, $xbs, lambda, t0, t1, $dispersions, $sources, $boundaries, nothing)
						end
					end
					continuousreleaseargs[2] = Symbol("t")
					qc.args[2].args[1].args = [qc.args[2].args[1].args[1]; continuousreleaseargs[1:end]...] # give it the correct set of arguments
					eval(qc)
					continuousreleaseargs[2] = Symbol("tau")
					qcf = quote
						function $(Symbol(string("long_", shortfunctionname, "_cf")))(thiswillbereplaced) # this function defines the continuous release function
							return kernel_cf(x, t, $x0s, $sigma0s, $vs, $sigmas, $Hs, $xbs, lambda, t0, t1, sourcestrength, $dispersions, $sources, $boundaries, nothing)
						end
					end
					continuousreleaseargs[2] = Symbol("t")
					qcf.args[2].args[1].args = [qcf.args[2].args[1].args[1]; continuousreleaseargs[1:end]...; :(sourcestrength::Function)] # give it the correct set of arguments
					eval(qcf)
					if documentation
						qdoc = quote
							@doc """$(DocumentFunction.documentfunction($(eval(Symbol(string("long_", shortfunctionname, "_ckernel")))), argtext=dictarguments, maintext="$($(numberofdimensions))-dimensional contaminant source kernel\n\n$($(docdispersions))\n\n$($(docsources))\n\n$($(docboundaries))"))""" $(Symbol(string("long_", shortfunctionname, "_ckernel")))
							@doc """$(DocumentFunction.documentfunction($(eval(Symbol(string("long_", shortfunctionname, "_c")))), argtext=dictarguments, maintext="$($(numberofdimensions))-dimensional continuous contaminant release with a unit flux\n\n$($(docdispersions))\n\n$($(docsources))\n\n$($(docboundaries))"))""" $(Symbol(string("long_", shortfunctionname, "_c")))
							@doc """$(DocumentFunction.documentfunction($(eval(Symbol(string("long_", shortfunctionname, "_cf")))), argtext=dictarguments, maintext="$($(numberofdimensions))-dimensional continuous contaminant release with a user-provided flux function\n\n$($(docdispersions))\n\n$($(docsources))\n\n$($(docboundaries))"))""" $(Symbol(string("long_", shortfunctionname, "_cf")))
							# @doc """$($(numberofdimensions))-dimensional contaminant source kernel \n- $($(docsources))\n- $($(docdispersions))\n- $($(docboundaries))\n$($(docarguments))""" $(Symbol(string("long_", shortfunctionname, "_ckernel")))
							# @doc """$($(numberofdimensions))-dimensional continuous contaminant release with a unit flux\n- $($(docsources))\n- $($(docdispersions))\n- $($(docboundaries))\n$($(docarguments))""" $(Symbol(string("long_", shortfunctionname, "_c")))
							# @doc """$($(numberofdimensions))-dimensional continuous contaminant release with a given flux\n- $($(docsources))\n- $($(docdispersions))\n- $($(docboundaries))\n$($(docarguments))""" $(Symbol(string("long_", shortfunctionname, "_cf")))
						end
						eval(qdoc)
					end
				end
			end
		end
	end
	MetaProgTools.replacesymbol!(bigq, :numberofdimensions, n)
	eval(bigq)
end

"""
Make documentation

$(DocumentFunction.documentfunction(documentationon))
"""
function documentationon()
	global documentation = true;
end

"""
Do not make documentation

$(DocumentFunction.documentfunction(documentationoff))
"""
function documentationoff()
	global documentation = false;
end

end
