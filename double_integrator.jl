using LinearAlgebra
using StaticArrays
using PyPlot
const Vec2f = SVector{2, Float64}
const Vec4f = SVector{4, Float64}
const s2f = Vec2f
const s4f = Vec4f


@inline function bisection_newton(f, df, left::Float64, right::Float64, eps=0.05, itr_max=20)::Float64
    x_est = (left + right)*0.5
    for itr = 1:itr_max
        f_df = f(x_est)/df(x_est)
        if (right-(x_est-f_df))*((x_est-f_df)-left)<0.0 || abs(f_df)<(right-left)/4.0
            f(x_est)>0 ? right = x_est : left = x_est
            f_df = x_est - (right + left)*0.5
        end
        x_est -= f_df
        abs(f_df)<eps && break
    end
    return x_est
end
"""
@inline function cost_optimal(s0::Vec4f, s1::Vec4f)
    a = Vec2f(s0[1:2])
    b = Vec2f(s0[3:4])
    c = Vec2f(s1[1:2])
    d =Vec2f(s1[3:4])
    return cost_optimal(a, b, c, d)
end
"""

@inline function cost_optimal(x0::Vec2f, v0::Vec2f, x1::Vec2f, v1::Vec2f)
    x01 = x1 - x0
    v01 = v1 - v0
    p = -4*(dot(v0, v0)+dot(v1, v1)+dot(v0, v1))
    q = 24*dot(v0+v1, x01)
    r = -36*dot(x01, x01)
    cost(t) = t + dot(v01, 4.0*v01/t-6(-v0*t+x01)/t^2)+dot(-6*v01/t^2+12(x01-v0*t)/t^3, -v0*t+x01)
    @fastmath d_cost(t) = t^4+p*t*t+q*t+r # df(t)/dt
    @fastmath dd_cost(t) = 4.0*t*t*t+2*p*t+q # ddf(t)/dt^2

    t_min = 0.0
    t_max = 10.0
    t_star = bisection_newton(d_cost, dd_cost, t_min, t_max)
    return cost(t_star)
end

# see ICRA paper: D.J.Webb et al Kinodynamic RRT* (2013)
function forward_reachable_box(x0::Vec2f, v0::Vec2f, r::Float64)
    @fastmath tau_x_plus = 2/3*(-(v0.^2).+r + v0.*sqrt.((v0.^2).+r))
    @fastmath tau_x_minus = 2/3*(-(v0.^2).+r - v0.*sqrt.((v0.^2).+r))
    @fastmath xmax = v0.*tau_x_plus + x0 + sqrt.(1/3*(tau_x_plus.^2).*(-tau_x_plus.+r))
    @fastmath xmin = v0.*tau_x_minus + x0 - sqrt.(1/3*(tau_x_minus.^2).*(-tau_x_minus.+r))

    @fastmath tau_v_plus = 0.5*r
    @fastmath vmax = v0 .+ sqrt.(tau_v_plus.*(-tau_v_plus.+r))
    @fastmath vmin = v0 .- sqrt.(tau_v_plus.*(-tau_v_plus.+r))
    return xmin, xmax, vmin, vmax
end

function backward_reachable_box(x0::Vec2f, v0::Vec2f, r::Float64)
    @fastmath tau_x_plus = 2/3*((v0.^2).-r + v0.*sqrt.((v0.^2).+r))
    @fastmath tau_x_minus = 2/3*((v0.^2).-r - v0.*sqrt.((v0.^2).+r))
    @fastmath xmax = v0.*tau_x_plus + x0 + sqrt.(1/3*(tau_x_plus.^2).*(tau_x_plus.+r))
    @fastmath xmin = v0.*tau_x_minus + x0 - sqrt.(1/3*(tau_x_minus.^2).*(tau_x_minus.+r))

    @fastmath tau_v_plus = 0.5*r
    @fastmath vmax = v0 .+ sqrt.(tau_v_plus.*(tau_v_plus.+r))
    @fastmath vmin = v0 .- sqrt.(tau_v_plus.*(tau_v_plus.+r))
    return xmin, xmax, vmin, vmax
end

function filter_freachable_exact(s_set, s_c, r)
    s_set_filtered = Vec4f[]
    for s in s_set
        @views xq = Vec2f(s[1:2])
        @views vq = Vec2f(s[3:4])
        @views xc = Vec2f(s_c[1:2])
        @views vc = Vec2f(s_c[3:4])
        cost_optimal(xc, vc, xq, vq)<r && push!(s_set_filtered, s)
    end
    return s_set_filtered
end

function filter_reachable(Sset::Vector{Vec4f}, idxset::Vector{Int64}, 
                          s_c::Vec4f, r::Float64, ForR::Symbol)
    s_set_filtered = Vec4f[]
    @views x_c = Vec2f(s_c[1:2])
    @views v_c = Vec2f(s_c[3:4])
    xmin, xmax, vmin, vmax = ((ForR==:F) ?
                              forward_reachable_box(x_c, v_c, r) :
                              backward_reachable_box(x_c, v_c, r))
                              
    function isinside(s::Vec4f)
        @inbounds !(xmin[1]<s[1]<xmax[1]) && return false
        @inbounds !(xmin[2]<s[2]<xmax[2]) && return false
        @inbounds !(vmin[1]<s[3]<vmax[1]) && return false
        @inbounds !(vmin[2]<s[4]<vmax[2]) && return false
        return true
    end

    idx_filtered = Int64[]
    for idx in idxset
        @inbounds s = Sset[idx]
        isinside(s) && push!(idx_filtered, idx)
        """
        if isinside(s)
            @views cost_optimal(s, s_c)<r && push!(idx_filtered, idx)
        end
        """
    end
    return idx_filtered
end

function test()
    N = 10^6
    s_set = Vec4f[]
    for i = 1:N
        push!(s_set, Vec4f(rand()-0.5, rand()-0.5, 3*rand()-1.5, 3*rand()-1.5))
    end
    s_c = Vec4f(+0, +0, 0.2, 0.2)
    idxset = [x for x in 1:N]
    idx_for = filter_reachable(s_set, idxset, s_c, 1.0, :F)
    idx_back = filter_reachable(s_set, idxset, s_c, 1.0, :B)
    s_for = s_set[idx_for]
    s_back = s_set[idx_back]

    setA = zeros(2, length(idx_for))
    setB = zeros(2, length(idx_back))
    for i = 1:length(idx_for)
        setA[:, i] = s_for[i][1:2]
    end
    for i = 1:length(idx_back)
        setB[:, i] = s_back[i][1:2]
    end
    scatter(setA[1, :], setA[2, :], s=3)
    scatter(setB[1, :], setB[2, :], s=3)
    xlim(-0.6, 0.6)
    ylim(-0.6, 0.6)
    return setA
end

