using LinearAlgebra
using StaticArrays
const Vec2f = SVector{2, Float64}
const Vec4f = SVector{4, Float64}

#using PyPlot
include("geometry.jl")
include("world.jl")
include("double_integrator.jl")


@inbounds @inline dist2(p, q)::Float64 = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)
# FMTree class
mutable struct FMTree
    s_init::Vec4f
    s_goal::Vec4f
    N #number of samples
    Pset::Vector{Vec4f}
    cost::Vector{Float64} #cost 
    time::Vector{Float64} #optimal time to connect one node to its root node
    parent::Vector{Int64} 
    bool_unvisit::BitVector #logical value for Vunvisit
    bool_open::BitVector #logical value for Open
    bool_closed::BitVector #logical value for Open
    world::World # simulation world config

    function FMTree(s_init::Vec4f, s_goal::Vec4f, N, world)
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
        time = zeros(N)
        parent = ones(Int, N)
        bool_unvisit = trues(N)
        bool_unvisit[1] = false
        bool_closed = falses(N)
        bool_open = falses(N)
        bool_open[1] = true
        new(s_init, s_goal,
            N, Pset, cost, time, parent, bool_unvisit, bool_open, bool_closed, world)
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
    idxset_tree = setdiff(union(idxset_open, idxset_closed), [1])
    """
    for idx in idxset_tree
        p1 = mat[:, idx]
        p2 = mat[:, this.parent[idx]]
        x = [p1[1], p2[1]]
        y = [p1[2], p2[2]]
        plot(x, y, c=:black, linewidth=1)
    end
    """
    for idx in idxset_tree
        s0 = this.Pset[this.parent[idx]]
        s1 = this.Pset[idx]
        tau = this.time[idx]
        show_trajectory(s0, s1, tau)
        println("fuck")
    end
    """
    scatter(mat[1, idxset_open], mat[2, idxset_open], c=:green, s=4)
    scatter(mat[1, idxset_closed], mat[2, idxset_closed], c=:black, s=5)
    scatter(mat[1, idxset_unvisit], mat[2, idxset_unvisit], c=:orange, s=5)
    xlim(this.world.x_min[1]-0.05, this.world.x_max[1]+0.05)
    ylim(this.world.x_min[2]-0.05, this.world.x_max[2]+0.05)
    """
end

function find_near_idx(Sset::Vector{Vec4f}, idxlst::Vector{Int64}, s_center::Vec4f, r::Float64)
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
    r = 1.0

    idxset_open = findall(this.bool_open)
    idxset_unvisit = findall(this.bool_unvisit)

    idx_lowest = idxset_open[findmin(this.cost[idxset_open])[2]]
    idxset_near, _, _ = filter_reachable(this.Pset, idxset_unvisit,
                                      this.Pset[idx_lowest], r, :F)

    for idx_near in idxset_near
        idxset_cand, distset_cand, timeset_cand = filter_reachable(this.Pset, idxset_open,
                                                     this.Pset[idx_near], r, :B)
        isempty(idxset_cand) && return
        cost_new, idx_costmin = findmin(this.cost[idxset_cand] + distset_cand)
        time_new = timeset_cand[idx_costmin] # optimal time for new connection
        idx_parent = idxset_cand[idx_costmin]
        if ~isIntersect(this.world, this.Pset[idx_near], this.Pset[idx_parent])
            this.bool_unvisit[idx_near] = false
            this.bool_open[idx_near] = true
            this.cost[idx_near] = cost_new
            this.time[idx_near] = time_new
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
t = FMTree(s_init, s_goal, 3000, wor)

@time for i=1:2000
    extend(t)
end








