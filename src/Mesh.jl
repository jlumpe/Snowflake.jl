module Mesh

using HexGrids


export make_mesh


# Define methods in MeshExt
function make_mesh end


function iter_tris(shape::HexagonShape)
	I = eltype(shape)
	a, b, c = hexaxes(I)
	shape2 = HexagonShape{I}(shape.n + 1)

	return Channel{NTuple{3, I}}() do ch
		for ix in shape2
			ix_a = ix + a
			ix_b = ix + b
			ix_ab = ix_a + b

			(ix_a in shape && ix_b in shape) || continue

			(ix in shape) && put!(ch, (ix, ix_a, ix_b))
			(ix_ab in shape) && put!(ch, (ix_b, ix_a, ix_ab))
		end
	end
end


default_cond = (ix) -> true


function linearize(shape::HexShape; cond=default_cond)
	lin2hex = eltype(shape)[]

	hex2lin = HexArray{Int}(shape)
	fill!(hex2lin, 0)

	next_i = 1
	for ix in shape
	   if cond(ix)
			push!(lin2hex, ix)
			hex2lin[ix] = next_i
			next_i += 1
		end
	end

	return (lin2hex, hex2lin)
end


end  # module Mesh
