using StaticArrays
using StaticArrays.ImmutableArrays
using PyPlot
export Polygon, isInside, isIntersect

const Vec2f = SVector{2, Float64}
const Vec4f = SVector{4, Float64}

abstract type Polygonic end

struct Polygon <: Polygonic
    N::Int 
    V::Vector{Vec2f}
    function Polygon(V_::AbstractVector) #V nust be counterclockwise
        N = length(V_)
        V = Vec2f[]
        for v in V_
            push!(V, Vec2f(v))
        end
        new(N, V)
    end
end

struct Rectangle <: Polygonic
    N::Int
    V::Vector{Vec2f}
    function Rectangle(center::T, width::Float64, height::Float64) where {T<:AbstractVector}
        V1 = Vec2f(center + [-width, -height]*0.5)
        V2 = Vec2f(center + [+width, -height]*0.5)
        V3 = Vec2f(center + [+width, +height]*0.5)
        V4 = Vec2f(center + [-width, +height]*0.5)
        new(4, [V1, V2, V3, V4])
    end
end


@inline function isInside(this::T, q::Vec2f) where {T<:Polygonic}
    # note: n'x = n'b
    q_vec = [q[1]; q[2]]
    for n = 1:this.N
        p1 = this.V[n]
        n<this.N ? p2 = this.V[n+1] : p2 = this.V[1]
        u = [p2[1]-p1[1]; p2[2]-p1[2]] #p2-p1
        n = [0 1; -1 0]*u # normla vector
        px = [q[1]-p1[1]; q[2]-p1[2]] # x-p1
        px'*n>0 && return false
    end
    return true
end


@inline function isIntersect(this::T, q1_::AbstractVector, q2_::AbstractVector) where {T<:Polygonic}
    @views q1 = Vec2f(q1_[1:2])
    @views q2 = Vec2f(q2_[1:2])
    return isIntersect(this, q1, q2)
end

@inline function isIntersect(this::T, q1::Vec2f, q2::Vec2f)::Bool where {T<:Polygonic}
    for n = 1:this.N
        p1 = this.V[n]
        if n<this.N
            p2=this.V[n+1]
        else
            p2=this.V[1]
        end
        isIntersect(p1, p2, q1, q2) && return true
    end
    return false
end

function show(this::T, color=:red) where {T<:Polygonic}
    for n = 1:this.N
        if n<this.N
            x = [this.V[n][1], this.V[n+1][1]]
            y = [this.V[n][2], this.V[n+1][2]]
        else
            x = [this.V[n][1], this.V[1][1]]
            y = [this.V[n][2], this.V[1][2]]
        end
        plot(x, y, "r-")
    end
end


@inline function isIntersect(p1::Vec2f, p2::Vec2f, q1::Vec2f, q2::Vec2f)::Bool
    # solve [p2-p1, q2-q1]*[s; t]=[q1-p1] or A*[s; t]=B
    # for computational efficiency, not using inv() func
    # A = [a, b; c, d]
    # B = [e; f]
    a = p2[1]-p1[1]
    b = -(q2[1]-q1[1])
    c = p2[2]-p1[2]
    d = -(q2[2]-q1[2])
    e = q1[1]-p1[1]
    f = q1[2]-p1[2]
    det = a*d-b*c
    abs(det)<1e-5 && return false

    s = 1/det*(d*e-b*f)
    t = 1/det*(-c*e+a*f)
    return (0.0<=s<=1.0) && (0.0<=t<=1.0)
end

function test()
    """
    v1 = [0.0, 0.0]
    v2 = [1.0, 0.0]
    v3 = [1.0, 1.0]
    v4 = [0.0, 1.0]
    V = [v1, v2, v3, v4]
    P = Polygon(V)
    """
    P = Rectangle([0.5 0.5], 1.0, 1.0)
    isInside(P, Vec2f([1.2, 0.8]))
end

