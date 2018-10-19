export Polygon, isInside, isIntersect

const Tuple2f = Tuple{Float64, Float64}

struct Polygon
    N::Int 
    V::Vector{Tuple2f}
    function Polygon(V) #V nust be counterclockwise
        N = length(V)
        new(N, V)
    end
end

@inline function isInside(this::Polygon, q::Tuple2f)
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

@inline function isIntersect(this::Polygon, q1::Tuple2f, q2::Tuple2f)::Bool
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
function show(this::Polygon, color=:red)
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


@inline function isIntersect(p1::Tuple2f, p2::Tuple2f, q1::Tuple2f, q2::Tuple2f)::Bool
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
    v1 = (0.0, 0.0)
    v2 = (1.0, 0.0)
    v3 = (1.0, 1.0)
    v4 = (0.0, 1.0)
    V = [v1, v2, v3, v4]
    P = Polygon(V)
    isInside(P, (1.2, 0.8))
end

