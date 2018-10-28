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

@inline function isIntersect(this::Sphere, x0::SVector{2, Float64}, x1::SVector{2, Float64})
    @fastmath v1 = this.c - x0
    @fastmath v2 = x1 - x0
    # calc dist between line segment and sphere orthogonal projection
    @fastmath vv = inpro(v1, v2)*v2/norm(v2)^2 # projected vectro
    @fastmath d = norm(x0 + vv - this.c)
    return d < this.r
end


""""
x_start = SVector(-2.0, 1.1)
x_end = SVector(2.0, 1.1)
s = Sphere(SVector(0., 0.), 1.2)
isCrossing(s, x_start, x_end)
"""
