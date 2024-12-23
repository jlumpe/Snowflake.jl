module MeshExt

using GeometryBasics

using HexGrids
using HexGrids: cartesian
using Snowflake: ReiterModel, frozen, boundary, nonreceptive
using Snowflake.Mesh: iter_tris, default_cond, linearize
import Snowflake.Mesh: make_mesh


default_zfunc = (ix) -> 0


function make_points(lin2hex::AbstractVector{<:HexIndex}; zfunc=default_zfunc, scale::Real=1)
	points = Point3f[]

	for (i, ix) in enumerate(lin2hex)
		x, y = cartesian(ix)
		z = convert(Float32, zfunc(ix))
		point = Point3f(x, y, z) * scale

		push!(points, point)
	end

	return points
end

make_points(lin2hex::AbstractVector{<:HexIndex}, a::HexArray) = make_points(lin2hex, zfunc=(ix) -> a[ix])


function make_tris(hex2lin::HexArray{<:Integer}, rev::Bool=false; offset::Int=0)
	tris = TriangleFace{Int}[]

	for (ix1, ix2, ix3) in iter_tris(keys(hex2lin))
		i1 = hex2lin[ix1]
		i2 = hex2lin[ix2]
		i3 = hex2lin[ix3]
		i1 > 0 && i2 > 0 && i3 > 0 || continue
		rev && ((i1, i2) = (i2, i1))
		push!(tris, TriangleFace{Int}(i1 + offset, i2 + offset, i3 + offset))
	end
	return tris
end


"""
Make point and triangle arrays.
"""
function make_mesh_parts(
		shape::HexShape;
		rev::Bool=false,
		zfunc=default_zfunc,
		cond=default_cond,
		scale::Real=1,
		tri_offset::Int=0,
	)
	lin2hex, hex2lin = linearize(shape, cond=cond)
	points = make_points(lin2hex, zfunc=zfunc, scale=scale)
	tris = make_tris(hex2lin, rev, offset=tri_offset);
	return (points, tris)
end

make_mesh_parts(array::HexArray; kw...) = make_mesh_parts(keys(array); kw..., zfunc=(ix) -> array[ix])


function make_mesh(
		model::ReiterModel;
		scale::Real=1.,
		zscale::Real=.1,
		zoffset::Real=0,
		normscale::Bool=true,
	)

	scale_::Float32 = scale
	zscale_::Float32 = zscale
	zoffset_::Float32 = zoffset

	if normscale
		scale_ /= model.shape.n
	else
		zscale_ *= model.shape.n
	end

	zarray = HexArray{Float32}(model.shape)
	for (ix, ct) in pairs(model.ct)
		if ct == frozen
			zarray[ix] = max(model.s[ix] + zoffset_, 0.f0) * zscale_
		else
			zarray[ix] = 0
		end
	end

	cond = (ix) -> (model.ct[ix] !== nonreceptive)

	bottom_points, bottom_tris = make_mesh_parts(
		model.shape,
		rev=true,
		cond=cond,
		scale=scale_,
	)
	top_points, top_tris = make_mesh_parts(
		zarray,
		cond=cond,
		scale=scale_,
		tri_offset=length(bottom_points),
	)

	# Concat
	points = vcat(bottom_points, top_points)
	tris = vcat(bottom_tris, top_tris)
	Mesh(points, tris)
end


end  # Module MeshExt
