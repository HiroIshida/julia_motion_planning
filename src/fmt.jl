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
    itr::Int64
    

    function FMTree(s_init::Vec4f, s_goal::Vec4f, N, world)
        println("initializing fmt ...")
        Pset = Vec4f[]
        push!(Pset, s_init) #inply idx_init = 1 
        myrn(min, max) = min + (max-min)*rand()
        for n = 1:N-2
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
        push!(Pset, s_goal) #inply idx_goal = N [last]
        cost = zeros(N)
        time = zeros(N)
        parent = ones(Int, N)
        bool_unvisit = trues(N)
        bool_unvisit[1] = false
        bool_closed = falses(N)
        bool_open = falses(N)
        bool_open[1] = true
        println("finish initializing")
        new(s_init, s_goal,
            N, Pset, cost, time, parent, bool_unvisit, bool_open, bool_closed, world, 0)
    end
end

function show(this::FMTree)
    println("drawing...")
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
    scatter(mat[1, 1], mat[2, 1], c=:blue, s=10, zorder = 100)
    scatter(mat[1, end], mat[2, end], c=:blue, s=10, zorder = 101)
    scatter(mat[1, idxset_open], mat[2, idxset_open], c=:gray, s=2)
    scatter(mat[1, idxset_closed], mat[2, idxset_closed], c=:gray, s=2)
    #scatter(mat[1, idxset_unvisit], mat[2, idxset_unvisit], c=:orange, s=5)
    for idx in idxset_tree
        s0 = this.Pset[this.parent[idx]]
        s1 = this.Pset[idx]
        tau = this.time[idx]
        show_trajectory(s0, s1, tau)
    end

    scatter(mat[1, 1], mat[2, 1], c=:blue, s=20, zorder = 100)
    scatter(mat[1, end], mat[2, end], c=:blue, s=20, zorder = 101)

    xlim(this.world.x_min[1]-0.05, this.world.x_max[1]+0.05)
    ylim(this.world.x_min[2]-0.05, this.world.x_max[2]+0.05)
    println("finish drawing")
end

function solve(this::FMTree, with_savefig = false)
    println("start solving")
    while(true)
        extend(this)
        if with_savefig
            if ((this.itr<100) & (this.itr % 20 == 1)) || (this.itr % 200==1)
                close()
                show(this)
                savefig("../fig/"*string(this.itr)*".png")
            end
        end
        !this.bool_unvisit[end] && break
    end
    idx = this.N
    idx_solution = Int64[idx]
    while(true)
        idx = this.parent[idx]
        push!(idx_solution, idx)
        idx == 1 && break
    end
    println("finish solving")
    return idx_solution
end

function extend(this::FMTree)
    this.itr += 1
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
        waypoints = gen_trajectory(this.Pset[idx_parent], this.Pset[idx_near], time_new, 10)
        if isValid(this.world, waypoints)
        #if ~isIntersect(this.world, this.Pset[idx_near], this.Pset[idx_parent])
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

