using StaticArrays
using StaticArrays.ImmutableArrays
const Vec2f = SVector{2, Float64}
const Vec4f = SVector{4, Float64}
@inline @inbounds norm(v) = sqrt(v[1]^2 + v[2]^2)
@inline @inbounds inpro(v1, v2) = v1[1]*v2[1]+v1[2]*v2[2]

abstract type GeoParametric end

struct Sphere <: GeoParametric
    c::SVector{2, Float64}
    r::Float64
end

@inline function isInside(this::Sphere, x::SVector{2, Float64})
    @fastmath d = norm(this.c-x)
    return (d < this.r )
end

@inline function isIntersect(this::Sphere, x_seq::Vector{Vec2f})
    for i in 1:length(x_seq)-1
        @inbounds x0 = x_seq[i]
        @inbounds x1 = x_seq[i+1]
        v1 = this.c - x0
        v2 = x1 - x0
        # calc dist between line segment and sphere orthogonal projection
        vv = inpro(v1, v2)*v2/norm(v2)^2 # projected vectro
        d = norm(x0 + vv - this.c)
        d < this.r && return true
    end
    return false
end



"""
x_0 = SVector(-1.0, 1.)
x_1 = SVector(0., 1.5)
x_2 = SVector(1.0, 1.)


s = Sphere(SVector(0., 0.), 1.3)
isIntersect(s, [x_0, x_1, x_2])
"""
