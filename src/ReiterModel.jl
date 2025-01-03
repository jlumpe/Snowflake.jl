@enum CellType frozen boundary nonreceptive

"""
α ∈ [0, 1) : Diffusion constant (or valid up to α=2?).
β ∈ [0, 1) : Amount of vater vapor in environment.
γ ≥ 0 : Constant amount added to receptive cells each step?
"""
const ReiterParams = @NamedTuple{α::Float64, β::Float64, γ::Float64}

const DEFAULT_PARAMS = (α=0, β=0, γ=0)


mutable struct ReiterModel
	const shape::HexagonShape{AxialIndex}
	params::ReiterParams
	"""Time index"""
	t::Int
	"""Cell states."""
	s::HexagonArray{Float64, Matrix{Float64}}
	"""Temporary space for calculation of next state array."""
	snext::HexagonArray{Float64, Matrix{Float64}}
	"""Cell types."""
	ct::HexagonArray{CellType, Matrix{CellType}}
	"""Side len of hex enclosing all "receptive" cells."""
	radius::Int

	function ReiterModel(n::Integer, params::ReiterParams=DEFAULT_PARAMS)
		model = new(
			HexagonShape(n),
			params,
			0,
			HexagonArray{Float64}(n),
			HexagonArray{Float64}(n),
			HexagonArray{CellType}(n),
		)
		init!(model)
	end
end

ReiterModel(n::Integer, params) = ReiterModel(n, ReiterParams(params))


"""
(Re)initialize the model state.

# Arguments
- `model`
- `n::Integer`: Size of initial hexagon of frozen cells at center.
- `random::Real`: If >0, init `s` of boundary cells with random value between 0 and this number.
- `rng`: RNG to use.
"""
function init!(
	model::ReiterModel,
	n::Integer=1;
	random::Real=0,
	rng::Random.AbstractRNG=Random.GLOBAL_RNG,
	)

	n < 0 && error("n must be nonnegative")
	randmax = Float32(random)

	model.t = 0
	model.radius = n

	fill!(model.s, model.params.β)
	fill!(model.ct, nonreceptive)

	for ix in HexagonShape(n + 1)
		ix in model.shape || continue
		if hexdist(ix) < n
			# Inner cell
			model.s[ix] = 1
			model.ct[ix] = frozen
		else
			# Boundary cell
			model.ct[ix] = boundary
			if random > 0
				model.s[ix] = rand(rng, Float32) * randmax
			end
		end
	end

	return model
end


"""
Refresh cell type array and `radius` attribute.

Not normally needed if starting from init! and calling update!.
"""
function refresh_ct!(model::ReiterModel)
	fill!(model.ct, nonreceptive)
	radius = 0

	for (z, s) in pairs(model.s)
		if s >= 1
			model.ct[z] = frozen
			radius = max(radius, hexdist(z))

			for z2 in neighbors(model.shape, z)
				if model.ct[z2] == nonreceptive
					model.ct[z2] = boundary
					radius = max(radius, hexdist(z2))
				end
			end
		end
	end

	model.radius = radius
	return model
end


"""
Update model by a single step.
"""
function update!(model::ReiterModel)
	# Update s values
	Threads.@threads for z in collect(model.shape)
		update_cell_s!(model, z)
	end
	(model.s, model.snext) = (model.snext, model.s)

	# Update cell states
	for z in model.shape
		model.ct[z] == nonreceptive || (model.radius = max(model.radius, hexdist(z)))

		# Boundary -> frozen transition
		if model.ct[z] === boundary && model.s[z] >= 1
			model.ct[z] = frozen
			# Update nonreceptive neighbors to boundary
			for zn in neighbors(model.shape, z)
				model.ct[zn] === nonreceptive && (model.ct[zn] = boundary)
			end
		end
	end

	model.t += 1

	return model
end


"""
Update s value of cell.
"""
function update_cell_s!(model::ReiterModel, z::AxialIndex)
	s = model.s[z]
	ct = model.ct[z]
	α, β, γ = model.params

	# Water that does (u) / does not (v) participate in diffusion
	u::Float64 = ct == nonreceptive ? s : 0.
	v::Float64 = ct == nonreceptive ? 0. : s

	# Diffiusion from neighboring cells
	neighboru = 0.
	for zn in neighbors(z)
		if zn in model.shape  # In grid
			# Only get diffusion from nonreceptive cells
			model.ct[zn] == nonreceptive && (neighboru += model.s[zn])
		else  # Past edge
			neighboru += β
		end
	end

	unext = u + α / 2 * (neighboru / 6 - u)
	vnext = ct == nonreceptive ? v : v + γ
	model.snext[z] = unext + vnext
end
