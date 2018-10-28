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
    """
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
    """
    for x in x_seq
        isInside(this, x) && return true
    end
    return false
end

function show(this::Sphere, color=:red)
    for theta in 0:0.3:2π
        x = this.c + this.r*Vec2f(cos(theta), sin(theta))
        scatter(x[1], x[2], c=color)
    end
end


"""
function test()
    x = Vec2f(-2., 1.6)
    x_seq = Vec2f[x]
    N = 30
    for i in 0:N
        dx = Vec2f(0.1, randn()*0.05)
        global x = x .+ dx
        push!(x_seq, x)
    end
    for i in 1:N
        x0 = x_seq[i]
        x1 = x_seq[i+1]
        plot([x0[1], x1[1]], [x0[2], x1[2]], c = :black)
    end


    x_0 = SVector(-1.0, 1.)
    x_1 = SVector(0., 1.5)
    x_2 = SVector(1.0, 1.)
    s = Sphere(SVector(0., 0.), 1.3)
    println(isIntersect(s, x_seq))
    show(s)
    println(s)
end
"""
