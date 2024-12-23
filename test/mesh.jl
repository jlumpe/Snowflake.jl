
using Snowflake.Mesh
using Snowflake.Mesh: linearize


function test_linearize(lin2hex, hex2lin)
	for (i, ix) in enumerate(lin2hex)
		@test hex2lin[ix] == i
	end

	for (ix, i) in pairs(hex2lin)
		@test i == 0 || lin2hex[i] == ix
	end
end


@testset "linearize" begin

	shape = HexagonShape(5)
	lin2hex, hex2lin = linearize(shape)
	test_linearize(lin2hex, hex2lin)

	cond(ix) = ix[1] % 2 == 0
	lin2hex, hex2lin = linearize(shape, cond=cond)
	test_linearize(lin2hex, hex2lin)

	@test length(lin2hex) > 0
	@test issetequal(lin2hex, collect(Iterators.filter(cond, shape)))
end
