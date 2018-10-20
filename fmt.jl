using LinearAlgebra
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
include("geometry.jl")
const Tuple2f = Tuple{Float64, Float64}
@inline dist2(p, q) = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)

function tupvec2mat(tupvec::Vector{Tuple2f})
    N = length(tupvec)
    mat = zeros(2, N)
    for i = 1:length(tupvec)
        for j=1:2
            mat[j, i] = tupvec[i][j]
        end
    end
    return mat
end


mutable struct World
    b_min
    b_max
    Pset::Vector{Polygon}
    function World(b_min, b_max, Pset)
        new(b_min, b_max, Pset)
    end
end
@inline function isValid(this::World, q::Tuple2f)
    ~((this.b_min[1]<q[1]<this.b_max[1])*(this.b_min[2]<q[2]<this.b_max[2])) && return false
    for P in this.Pset
        isInside(P, q) && return false
    end
    return true
end
@inline function isIntersect(this::World, q1::Tuple2f, q2::Tuple2f)
    for P in this.Pset
        isIntersect(P, q1, q2) && return true
    end
    return false
end

function show(this::World)
    p1 = [this.b_min[1], this.b_min[2]]
    p2 = [this.b_min[1], this.b_max[2]]
    p3 = [this.b_max[1], this.b_max[2]]
    p4 = [this.b_max[1], this.b_min[2]]
    plot([p1[1], p2[1]], [p1[2], p2[2]], "r-")
    plot([p2[1], p3[1]], [p2[2], p3[2]], "r-")
    plot([p3[1], p4[1]], [p3[2], p4[2]], "r-")
    plot([p4[1], p1[1]], [p4[2], p1[2]], "r-")
    for P in this.Pset
        show(P)
    end
end



# FMTree class
mutable struct FMTree
    x_init::Tuple2f
    x_goal::Tuple2f
    N #number of samples
    Pset::Vector{Tuple2f}
    cost::Vector{Float64} #cost 
    parent::Vector{Int64} 
    bool_unvisited::BitVector #logical value for Vunvisited
    bool_open::BitVector #logical value for Open
    bool_closed::BitVector #logical value for Open
    world::World # simulation world config

    function FMTree(x_init, x_goal, N, world)
        Pset = Tuple2f[]
        push!(Pset, x_init)
        for n = 2:N
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
        new(x_init, x_goal,
            N, Pset, cost, parent, bool_unvisited, bool_open, bool_closed, world)
    end
end

function show(this::FMTree)
    show(this.world)
    mat =tupvec2mat(this.Pset)
    idxset_open = findall(this.bool_open)
    idxset_closed = findall(this.bool_closed)
    idxset_unvisited = findall(this.bool_unvisited)
    idxset_tree = union(idxset_open, idxset_closed)
    for idx in idxset_tree
        p1 = mat[:, idx]
        p2 = mat[:, this.parent[idx]]
        x = [p1[1], p2[1]]
        y = [p1[2], p2[2]]
        plot(x, y, c=:black, linewidth=1)
    end

    scatter(mat[1, idxset_open], mat[2, idxset_open], c=:green, s=4)
    scatter(mat[1, idxset_closed], mat[2, idxset_closed], c=:black, s=5)
    scatter(mat[1, idxset_unvisited], mat[2, idxset_unvisited], c=:orange, s=5)
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
    for idx in idxset_unvisited
        @inbounds dist = dist2(this.Pset[idx_lowest],this.Pset[idx])
        if dist<radi
            push!(idxset_neighbor, idx)
        end
    end
   
    #find 
    for idx_neighbor in idxset_neighbor
        #serach candidate points for connection
        idxset_candidate = Int64[];
        distset_candidate = Float64[];
        for idx_open in idxset_open
            @inbounds dist = dist2(this.Pset[idx_open], this.Pset[idx_neighbor])
            if dist<radi
                push!(idxset_candidate, idx_open)
                push!(distset_candidate, dist)
            end
        end
        isempty(idxset_candidate) && return

        #cost-optimal connection
        tmp = findmin(this.cost[idxset_candidate] + distset_candidate)
        cost_new = tmp[1];
        idx_parent = idxset_candidate[tmp[2]]
        if ~isIntersect(this.world, this.Pset[idx_neighbor], this.Pset[idx_parent])
            this.bool_unvisited[idx_neighbor] = false
            this.bool_open[idx_neighbor] = true
            this.cost[idx_neighbor] = cost_new
            this.parent[idx_neighbor] = idx_parent
        end
    end

    this.bool_open[idx_lowest] = false
    this.bool_closed[idx_lowest] = true
end

v1 = (0.25, 0.25)
v2 = (0.5, 0.5)
v3 = (0.2, 0.5)

Pset = [Polygon([v1, v2, v3])]
wor = World([0, 0], [1.0, 1.0], Pset)
t = FMTree((0.1, 0.1), (0.9, 0.9), 10000, wor)

@time for i=1:10000
    println(i)
    extend(t)
end









