using LinearAlgebra
using StaticArrays
using StaticArrays.ImmutableArrays
#using PyPlot
include("geometry.jl")
const Vec2f = SVector{2, Float64}
const Vec4f = SVector{4, Float64}

@inbounds @inline dist2(p, q)::Float64 = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)

mutable struct World
    x_min
    x_max
    v_min
    v_max
    Pset::Vector{Polygonic}
    function World(x_min, x_max, v_min, v_max, Pset)
        new(x_min, x_max, v_min, v_max, Pset)
    end
end

@inline function isValid(this::World, s_q::Vec4f)
    !(this.x_min[1]<s_q[1]<this.x_max[1]) && return false
    !(this.x_min[2]<s_q[2]<this.x_max[2]) && return false
    !(this.v_min[1]<s_q[3]<this.v_max[1]) && return false
    !(this.v_min[2]<s_q[4]<this.v_max[2]) && return false
    for P in this.Pset
        isInside(P, Vec2f(s_q[1:2])) && return false
    end
    return true
end

@inline function isIntersect(this::World, q1::Vec4f, q2::Vec4f)
    for P in this.Pset
        isIntersect(P, q1, q2) && return true
    end
    return false
end

function show(this::World)
    p1 = [this.x_min[1], this.x_min[2]]
    p2 = [this.x_min[1], this.x_max[2]]
    p3 = [this.x_max[1], this.x_max[2]]
    p4 = [this.x_max[1], this.x_min[2]]
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
    s_init::Vec4f
    s_goal::Vec4f
    N #number of samples
    Pset::Vector{Vec4f}
    cost::Vector{Float64} #cost 
    parent::Vector{Int64} 
    bool_unvisit::BitVector #logical value for Vunvisit
    bool_open::BitVector #logical value for Open
    bool_closed::BitVector #logical value for Open
    world::World # simulation world config
    metric

    function FMTree(s_init::Vec4f, s_goal::Vec4f, N, world, metric)
        Pset = Vec4f[]
        push!(Pset, s_init)
        myrn(min, max) = min + (max-min)*rand()
        for n = 2:N
            while(true)
                p = Vec4f(myrn(world.x_min[1], world.x_max[1]), 
                          myrn(world.x_min[2], world.x_max[2]),
                          myrn(world.v_min[1], world.v_max[1]),
                          myrn(world.v_min[1], world.v_max[1]))
                if isValid(world, p)
                    push!(Pset, p)
                    break
                end
            end
        end
        cost = zeros(N)
        parent = ones(Int, N)
        bool_unvisit = trues(N)
        bool_closed = falses(N)
        bool_open = falses(N)
        bool_open[1] = true
        new(s_init, s_goal,
            N, Pset, cost, parent, bool_unvisit, bool_open, bool_closed, world, metric)
    end
end

function show(this::FMTree)
    show(this.world)
    N = length(this.Pset)
    mat = zeros(2, N)
    for idx = 1:N
        mat[:, idx] = this.Pset[idx][1:2]
    end
    idxset_open = findall(this.bool_open)
    idxset_closed = findall(this.bool_closed)
    idxset_unvisit = findall(this.bool_unvisit)
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
    scatter(mat[1, idxset_unvisit], mat[2, idxset_unvisit], c=:orange, s=5)
    xlim(this.world.x_min[1]-0.05, this.world.x_max[1]+0.05)
    ylim(this.world.x_min[2]-0.05, this.world.x_max[2]+0.05)
end

function cleanup(this::FMTree)
    this.cost = zeros(this.N)
    this.bool_open = falses(this.N); this.bool_open[1] = true
    this.bool_unvisit = trues(this.N)
    this.bool_closed = falses(this.N)
end

function find_near_idx(s_center::Vec4f, Sset::Vector{Vec4f}, idxlst::Vector{Int64}, r::Float64)
    idxset_near = Int64[]
    distset_near = Float64[]
    for idx in idxlst
        @inbounds dist = dist2(s_center, Sset[idx])
        if dist<r
            push!(idxset_near, idx)
            push!(distset_near, dist)
        end
    end
    return idxset_near, distset_near
end

function find_near_idx(s_center::Vec4f, Sset::Vector{Vec4f}, idxlst::Vector{Int64}, r::Float64, ForR::Symbol)
    idxset_near = Int64[]
    distset_near = Float64[]
    for idx in idxlst
        @inbounds dist = dist2(s_center, Sset[idx])
        if dist<r
            push!(idxset_near, idx)
            push!(distset_near, dist)
        end
    end
    return idxset_near, distset_near
end


function extend(this::FMTree)
    r = 0.1

    #select the lowest-cost node 
    idxset_open = findall(this.bool_open)
    idx_lowest = idxset_open[findmin(this.cost[idxset_open])[2]]

    #find nears within Vunvisit
    idxset_unvisit = findall(this.bool_unvisit)
    idxset_near, distset_near = find_near_idx(this.Pset[idx_lowest], this.Pset, idxset_unvisit, r)
   
    #find 
    for idx_near in idxset_near
        #serach cand points for connection
        idxset_cand, distset_cand = find_near_idx(this.Pset[idx_near], this.Pset, idxset_open, r)
        isempty(idxset_cand) && return

        #cost-optimal connection
        tmp = findmin(this.cost[idxset_cand] + distset_cand)
        cost_new = tmp[1];
        idx_parent = idxset_cand[tmp[2]]
        if ~isIntersect(this.world, this.Pset[idx_near], this.Pset[idx_parent])
            this.bool_unvisit[idx_near] = false
            this.bool_open[idx_near] = true
            this.cost[idx_near] = cost_new
            this.parent[idx_near] = idx_parent
        end
    end

    this.bool_open[idx_lowest] = false
    this.bool_closed[idx_lowest] = true
end



v1 = (0.25, 0.25)
v2 = (0.5, 0.5)
v2 = (0.2, 0.5)

#Pset = [Polygon([[0.2, 0.2], [0.4, 0.2], [0.3, 0.3]])]
Pset = [Rectangle([0.2, 0.5], 0.05, 0.8),
        Rectangle([0.4, 0.5], 0.05, 0.8)]
        
x_min = [0, 0]
x_max = [1.0, 1.0]
v_min = [-0.5, -0.5]
v_max = [0.5, 0.5]
wor = World(x_min, x_max, v_min, v_max, Pset)

s_init = Vec4f([0.1, 0.1, 0.0, 0.0])
s_goal = Vec4f([0.9, 0.9, 0.0, 0.0])
t = FMTree(s_init, s_goal, 3000, wor, dist2)

@time for i=1:2000
    extend(t)
end








