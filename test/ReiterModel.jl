using Snowflake: ReiterModel, frozen, boundary, nonreceptive
using Snowflake: init!, refresh_ct!


@testset "init" begin
	params = (α=2, β=.4, γ=.001)
	model = ReiterModel(30, params)

	n = 3
	init!(model, n)

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
