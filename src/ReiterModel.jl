@enum CellType frozen boundary nonreceptive

"""
α ∈ [0, 1) : Diffiusion constant (or valid up to α=2?).
β ∈ [0, 1) : Amount of vater vapor in environment.
γ ≥ 0 : Constant amount added to receptive cells each step?
"""
const ReiterParams = @NamedTuple{α::Float64, β::Float64, γ::Float64}


mutable struct ReiterModel
	const shape::HexagonShape{AxialIndex}
	params::ReiterParams
	t::Int
	s::HexagonHexArray{Float64, AxialIndex}
	snext::HexagonHexArray{Float64, AxialIndex}
	ct::HexagonHexArray{CellType, AxialIndex}

	function ReiterModel(n::Int, params=(α=0, β=0, γ=0))
		model = new(HexagonShape(n), params, 0, HexagonHexArray{Float64}(n), HexagonHexArray{Float64}(n), HexagonHexArray{CellType}(n))
		init!(model)
	end
end


"""
(Re)initialize the model state.
"""
function init!(model::ReiterModel)
	model.t = 0

	o = AxialIndex()
	fill!(model.s, model.params.β)
	model.s[o] = 1

	fill!(model.ct, nonreceptive)
	model.ct[o] = frozen
	for z in neighbors(o)
		model.ct[z] = boundary
	end

	return model
end


"""
Refresh cell type array.

Not normally needed if starting from init! and calling update!.
"""
function refresh_ct!(model::ReiterModel)
	frozencells = AxialIndex[]

	for (z, s) in pairs(model.s)
		if s >= 1
			model.ct[z] = frozen
			push!(frozencells, z)
		else
			model.ct[z] = nonreceptive
		end
	end

	for z in frozencells
		for nz in neighbors(model.shape, z)
			model.ct[nz] === nonreceptive && (model.ct[nz] = boundary)
		end
	end

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
