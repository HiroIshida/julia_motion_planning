using LinearAlgebra
using StaticArrays
using PyPlot
const Vec2f = SVector{2, Float64}
const Vec4f = SVector{4, Float64}

include("src/geometry.jl")
include("src/world.jl")
include("src/double_integrator.jl")
include("src/fmt.jl")


obstacle_set = [Rectangle([0.2, 0.3], 0.2, 0.2),
        Rectangle([0.5, 0.5], 0.2, 0.3),
        Rectangle([0.8, 0.3], 0.2, 0.1), 
        Rectangle([0.8, 0.6], 0.15, 0.2),
        Rectangle([0.2, 0.7], 0.1, 0.4),
        Rectangle([0.2, 0.7], 0.1, 0.4)]
x_min = [0, 0]
x_max = [1.0, 1.0]
v_min = [-0.5, -0.5]
v_max = [0.5, 0.5]
W = World(x_min, x_max, v_min, v_max, obstacle_set)

s_init = Vec4f([0.1, 0.1, 0.0, 0.0])
s_goal = Vec4f([0.9, 0.9, 0.0, 0.0])

Nsample = 3000
fmt = @time FMTree(s_init, s_goal, Nsample, W)

with_savefig = false
idx_solution = @time solve(fmt, with_savefig)

show(fmt)
for idx in idx_solution
    s0 = fmt.Pset[fmt.parent[idx]]
    s1 = fmt.Pset[idx]
    tau = fmt.time[idx]
    show_trajectory(s0, s1, tau, 20, :blue, 1.5)
end
savefig("./fig/final.png")



