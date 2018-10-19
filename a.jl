@inline dist2(p, q) = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)

# FMTree class
mutable struct FMTree
    N #number of samples
    Pset #points
    bool_unvisited #logical value for Vunvisited
    function FMTree(N)
        Pset = Array{Float64, 2}[]
        for n = 1:N
            p_cand = rand(2, 1);
            push!(Pset, p_cand)
        end
        bool_unvisited = trues(N)
        new(N, Pset, bool_unvisited)
    end
end

function extend(this::FMTree)
    #find neighbors within Vunvisited
    idxset_unvisited = findall(this.bool_unvisited)
    bool = this.bool_unvisited;
    idxset_unvisited = findall(bool)
    for idx in idxset_unvisited
        dist = dist2([0 0], this.Pset[idx])
    end
end

t = FMTree(10000)
@time for i=1:10000
    println(i)
    extend(t)
end
