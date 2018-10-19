using LinearAlgebra
const Tuple2f = Tuple{Float64, Float64}
@inline dist2(p, q) = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)
# World 
mutable struct World
    b_min
    b_max
    function World(b_min, b_max)
        new(b_min, b_max)
    end
end
@inline function isValid(this::World, p)
    return (this.b_min[1]<p[1]<this.b_max[1])
    *(this.b_min[2]<p[2]<this.b_max[2])
end

# FMTree class
mutable struct FMTree
    N #number of samples
    Pset::Vector{Tuple2f}
    cost::Vector{Float64} #cost 
    parent::Vector{Int64} 
    bool_unvisited::BitVector #logical value for Vunvisited
    bool_open::BitVector #logical value for Open
    bool_closed::BitVector #logical value for Open
    world::World # simulation world config

    function FMTree(N, world)
        Pset = Tuple2f[]
        for n = 1:N
            while(true)
                p = (rand(), rand())
                if isValid(world, p)
                    push!(Pset, p)
                    break
                end
            end
        end
        cost = zeros(N)
        parent = ones(Int, N)
        bool_unvisited = trues(N)
        bool_closed = falses(N)
        bool_open = falses(N)
        bool_open[1] = true
        new(N, Pset, cost, parent, bool_unvisited, bool_open, bool_closed, world)
    end
end

function show(this::FMTree)
    idxset_open = findall(this.bool_open)
    idxset_closed = findall(this.bool_closed)
    idxset_unvisited = findall(this.bool_unvisited)
    idxset_tree = union(idxset_open, idxset_closed)
    for idx in idxset_tree
        p1 = this.Pset[idx]
        p2 = this.Pset[this.parent[idx]]
        x = [p1[1], p2[1]]
        y = [p1[2], p2[2]]
        plot(x, y, c=:black)
    end

    for i in idxset_closed
        scatter(this.Pset[i][1], this.Pset[i][2], c=:orange)
    end
    for i in idxset_open
        scatter(this.Pset[i][1], this.Pset[i][2], c=:green)
    end
    for i in idxset_unvisited
        scatter(this.Pset[i][1], this.Pset[i][2], c=:red)
    end
    xlim(-0.1, 1.1)
    ylim(-0.1, 1.1)
end

function cleanup(this::FMTree)
    this.cost = zeros(this.N)
    this.bool_open = falses(this.N); this.bool_open[1] = true
    this.bool_unvisited = trues(this.N)
    this.bool_closed = falses(this.N)
end


function extend(this::FMTree)
    radi = 0.05

    #select the lowest-cost node 
    idxset_open = findall(this.bool_open)
    idx_lowest = idxset_open[findmin(this.cost[idxset_open])[2]]

    #find neighbors within Vunvisited
    idxset_unvisited = findall(this.bool_unvisited)
    idxset_neighbor = Int64[]
    s=0
    for idx in idxset_unvisited
        @inbounds dist = dist2(this.Pset[idx_lowest],this.Pset[idx])
        if dist<radi
            push!(idxset_neighbor, idx)
        end
        s += 1
    end
   
    #find 
    for idx_neighbor in idxset_neighbor
        #serach candidate points for connection
        idxset_candidate = [];
        distset_candidate = [];
        for idx_open in idxset_open
            @inbounds dist = dist2(this.Pset[idx_open], this.Pset[idx_neighbor])
            if dist<radi
                push!(idxset_candidate, idx_open)
                push!(distset_candidate, dist)
            end
        end

        #cost-optimal connection
        if ~isempty(idxset_candidate)
            tmp = findmin(this.cost[idxset_candidate] + distset_candidate)
            cost_new = tmp[1];
            idx_parent = idxset_candidate[tmp[2]]
            this.bool_unvisited[idx_neighbor] = false
            this.bool_open[idx_neighbor] = true
            this.cost[idx_neighbor] = cost_new
            this.parent[idx_neighbor] = idx_parent
            
        end
    end
    this.bool_open[idx_lowest] = false
    this.bool_closed[idx_lowest] = true
end

wor = World([0, 0], [1.0, 1.0])
t = FMTree(10000, wor)

@time for i=1:10000
    println(i)
    extend(t)
end








