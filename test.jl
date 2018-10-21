@inline dist2(p, q) = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)
function rand_gen()
    r2set = Tuple{Float64,Float64}[]
    for i=1:10000
        r2_add = (rand(), rand())
        push!(r2set, r2_add)
    end
    return r2set
end

using StaticArrays
struct mytype
    Vector{SVector}[]
end

function test()
    N = 10000
    r2set = rand_gen()
    a = (1,1)
    b = (2,2)

    s = 0.0
    @time for i=1:N, j=1:N
        @inbounds s += dist2(r2set[i], r2set[j])
    end

    @time for i=1:N, j=1:N
        s += dist2(a, b)
    end
end
test()
