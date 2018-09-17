import QuadGK

include("gencode.jl")

function kernel(x::Vector, tau::Real, x0::Vector, sigma0::Vector, v::Vector, sigma::Vector, H::Vector, xb::Vector, dispersions, sources, boundaries, distributions=nothing)
	return innerkernel(Val{length(x)}, x, tau, x0, sigma0, v, sigma, H, xb, dispersions, sources, boundaries, distributions)
end

function cinnerkernel(dimtype, x, tau, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, t, dispersions, sources, boundaries, distributions)
	if inclosedinterval(t - tau, t0, t1)
		return exp(-lambda * tau) * innerkernel(dimtype, x, tau, x0, sigma0, v, sigma, H, xb, dispersions, sources, boundaries, distributions)
	else
		return 0.
	end
end

function quadgkwithtol(f, a, b)
	return QuadGK.quadgk(f, a, b; rtol=1.0e-7, atol=1e-8)[1]
end

function kernel_c(x::Vector, t::Real, x0::Vector, sigma0::Vector, v::Vector, sigma::Vector, H::Vector, xb::Vector, lambda::Real, t0::Real, t1::Real, dispersions, sources, boundaries, distributions=nothing)
	dimtype = Val{length(x)}
	return kernel_c(dimtype, x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, dispersions, sources, boundaries, distributions)
end
function kernel_c(dimtype, x::Vector, t::Real, x0::Vector, sigma0::Vector, v::Vector, sigma::Vector, H::Vector, xb::Vector, lambda::Real, t0::Real, t1::Real, dispersions, sources, boundaries, distributions=nothing)
	return kernel_cf(dimtype, x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, t->inclosedinterval(t, t0, t1) ? 1. : 0., dispersions, sources, boundaries, distributions)
end

function kernel_cf(x::Vector, t::Real, x0::Vector, sigma0::Vector, v::Vector, sigma::Vector, H::Vector, xb::Vector, lambda::Real, t0::Real, t1::Real, sourcestrength::Function, dispersions, sources, boundaries, distributions=nothing)
	dimtype = Val{length(x)}
	return kernel_cf(dimtype, x, t, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, sourcestrength, dispersions, sources, boundaries, distributions)
end
function kernel_cf(dimtype, x::Vector, t::Real, x0::Vector, sigma0::Vector, v::Vector, sigma::Vector, H::Vector, xb::Vector, lambda::Real, t0::Real, t1::Real, sourcestrength::Function, dispersions, sources, boundaries, distributions=nothing)
	if t - t0 <= 0
		return 0.0
	elseif t - t1 <= 0 && inclosedinterval(t - t0, 0, t)
		return quadgkwithtol(tau->sourcestrength(t - tau) * cinnerkernel(dimtype, x, tau, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, t, dispersions, sources, boundaries, distributions), 0, t - t0)
	elseif 0 <= t - t1 && t - t0 <= t
		return quadgkwithtol(tau->sourcestrength(t - tau) * cinnerkernel(dimtype, x, tau, x0, sigma0, v, sigma, H, xb, lambda, t0, t1, t, dispersions, sources, boundaries, distributions), t - t1, t - t0)
	elseif inclosedinterval(t - t1,0,t) && t - t0 >= t
		error("t0 is less than zero, but the code assumes that t0>=0")
	else
		error("outside of ifelses: [t, t0, t1] = [$(t), $(t0), $(t1)]")
	end
end

@gen_code function innerkernel(dimtype::Type{Val{dimensions}}, x::Vector, tau::Real, x0::Vector, sigma0::Vector, v::Vector, sigma::Vector, H::Vector, xb::Vector, ::Type{Val{dispersions}}, ::Type{Val{sources}}, ::Type{Val{boundaries}}, distributions) where {dimensions,dispersions,sources,boundaries}
	if dimensions <= 0
		error("Dimensions must be positive")
	end
	distexprs = distributionexprs(distributions, dimensions)
	for i = 1:dimensions
		@code distexprs[i]
		@code dispersiontimedependenceexpr(dispersions[i], sources[i], i)
		ide = infinitedomainexpr(sources[i], i)
		be = boundaryexpr(ide, boundaries[i], i)
		@code :($(symbolindex("retval", i)) = $be)
	end
	@code retexpr(dimensions)
	return code
end

function retexpr(dimensions)
	q = :(retval1)
	for i = 2:dimensions
		q = :($q * $(symbolindex("retval", i)))
	end
	return :(return $q)
end

function distributionexprs(distributions::Type{Nothing}, dimensions::Int)
	distexprs = Array{Expr}(undef, dimensions)
	for i = 1:dimensions
		distexprs[i] = :($(Symbol(string("dist", i))) = Anasol.standardnormal)
	end
	return distexprs
end

function distributionexprs(distributions, dimensions)
	distexprs = Array{Expr}(undef, dimensions)
	for i = 1:dimensions
		distexprs[i] = :($(Symbol(string("dist", i))) = distributions[$i])
	end
	return distexprs
end

function symbolindex(s, i::Int)
	return Symbol(string(s, i))
end

function getinitdispersionfactor(sourcetype, i::Int)
	if sourcetype == :box
		return :(0)
	elseif sourcetype == :dispersed
		return :(sigma0[$i] ^ 2)
	else
		error("Unknown source type: $sourcetype")
	end
end

function dispersiontimedependenceexpr(D, sourcetype, i::Int)
	initdispersionexp = getinitdispersionfactor(sourcetype, i)
	if D == :linear
		timedispersionexp = :(sigma[$i] ^ 2 * tau)
		dispersionexp = :(sqrt($initdispersionexp + $timedispersionexp))
	elseif D == :fractional
		timedispersionexp = :(sigma[$i] ^ 2 * tau ^ (2 * H[$i]))
		dispersionexp = :(sqrt($initdispersionexp + $timedispersionexp))
	elseif typeof(D) == Expr
		dispersionexp = D
	else
		error("Unknown dispersion time dependence $D")
	end
	return :($(symbolindex("sigmat", i)) = $dispersionexp)
end

function infinitedomainexpr(S, i)
	distsym = symbolindex("dist", i)
	sigmatsym = symbolindex("sigmat", i)
	if S == :dispersed
		return :(Distributions.pdf($distsym,((point - x0[$i]) - v[$i] * tau) / $sigmatsym) / $sigmatsym)
	elseif S == :box
		return :((Distributions.cdf($distsym,(((point - x0[$i]) - v[$i] * tau) + 0.5 * sigma0[$i]) / $sigmatsym) - Distributions.cdf($distsym,(((point - x0[$i]) - v[$i] * tau) - 0.5 * sigma0[$i]) / $sigmatsym)) / sigma0[$i])
	else
		error("Unknown source type: $S")
	end
end

function boundaryexpr(infdomexpr, B, i)
	if B == :infinite
		return MetaProgTools.replacesymbol(infdomexpr, :point, :(x[$i]))
	elseif B == :reflecting || B == :absorbing
		e1 = MetaProgTools.replacesymbol(infdomexpr, :point, :(x[$i]))
		e2 = MetaProgTools.replacesymbol(infdomexpr, :point, :(2 * xb[$i] - x[$i]))
		if B == :reflecting
			return :($e1 + $e2)
		else
			return :($e1 - $e2)
		end
	else
		error("Unknown boundary condition $B")
	end
end
