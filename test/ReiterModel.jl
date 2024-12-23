using Snowflake: ReiterModel, frozen, boundary, nonreceptive
using Snowflake: init!, update!, refresh_ct!


@testset "init" begin
	params = (α=2, β=.4, γ=.001)
	model = ReiterModel(30, params)

	n = 3
	init!(model, n)
	@test model.radius == n

	for ix in model.shape
		d = hexdist(ix)
		if d < n
			@test model.s[ix] == 1
			@test model.ct[ix] == frozen
		elseif d == n
			@test model.s[ix] == params.β
			@test model.ct[ix] == boundary
		else
			@test model.s[ix] == params.β
			@test model.ct[ix] == nonreceptive
		end
	end
end


@testset "run" begin
	params = (α=2, β=.4, γ=.001)
	model = ReiterModel(30, params)
	init!(model)

	for _ in 1:50
		update!(model)
	end

	radius = 0

	# Check cells
	for z in model.shape
		ct = model.ct[z]
		s = model.s[z]

		# Check frozen based on value of s
		@test (ct == frozen) == (s >= 1)

		ct == nonreceptive || (radius = max(radius, hexdist(z)))

		# Check neighbors
		has_frozen_neighbor = false
		for zn in neighbors(model.shape, z)
			ctn = model.ct[zn]

			# Test neighbors of frozen are frozen or boundary
			if ct == frozen
				@test ctn in (frozen, boundary)
			end

			ctn == frozen && (has_frozen_neighbor = true)
		end

		# Boundary cells must have frozen neighbor
		if ct == boundary
			@test has_frozen_neighbor
		end
	end

	@test model.radius == radius
end


@testset "refresh_ct!()" begin
	params = (α=2, β=.4, γ=.001)
	model = ReiterModel(30, params)
	init!(model)

	for _ in 1:50
		update!(model)
	end

	old_ct = copy(model.ct)
	old_radius = model.radius

	fill!(model.ct, nonreceptive)
	model.radius = 0
	refresh_ct!(model)

	@test model.ct == old_ct
	@test model.radius == old_radius
end
