using LinearAlgebra
using PyPlot

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
    Pset #points
    cost #cost 
    parent 
    bool_unvisited #logical value for Vunvisited
    bool_open #logical value for Open
    bool_closed #logical value for Open
    world # simulation world config

    function FMTree(N, world)
        Pset = Array{Float64, 2}[]
        for n = 1:N
            while(true)
                p_cand = rand(2, 1);
                if isValid(world, p_cand)
                    push!(Pset, p_cand)
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
    idx_open = findall(this.bool_open)
    idx_closed = findall(this.bool_closed)
    idx_unvisited = findall(this.bool_unvisited)
    idxset_tree = union(idx_open, idx_closed)
    for idx in idxset_tree
        p1 = this.Pset[idx]
        p2 = this.Pset[this.parent[idx]]
        x = [p1[1], p2[1]]
        y = [p1[2], p2[2]]
        plot(x, y, c=:black)
    end

    scatter(this.Pset[idx_closed][1], this.Pset[idx_closed][2], c=:orange)
    scatter(this.Pset[idx_open][1], this.Pset[idx_open][2], c=:green)
    scatter(this.Pset[idx_unvisited][1], this.Pset[idx_unvisited][2], c=:red)
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
    bool = this.bool_unvisited;
    bool = trues(10000)
    idxset_unvisited = findall(bool)
    idxset_neighbor = Int64[]
    for idx in idxset_unvisited
        """
        dist = dist2(this.Pset[idx_lowest], this.Pset[idx])
        if dist<radi
            push!(idxset_neighbor, idx)
        end
        """
    end
   
    """
    #find 
    for idx_neighbor in idxset_neighbor
        #serach candidate points for connection
        idxset_candidate = [];
        distset_candidate = [];
        for idx_open in idxset_open
            dist = dist2(this.Pset[idx_open], this.Pset[idx_neighbor])
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
    """
end

wor = World([0, 0], [1.0, 1.0])
t = FMTree(10000, wor)

@time for i=1:1
    println(i)
    extend(t)
end








