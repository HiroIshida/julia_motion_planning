using LinearAlgebra
using StaticArrays
using PyPlot
const SVector2f = SVector{2, Float64}
const SVector4f = SVector{4, Float64}
const s2f = SVector2f
const s4f = SVector4f


function bisection_newton(f, df, left::Float64, right::Float64, eps=0.05, itr_max=20)::Float64
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

function find_tau_star(x0::SVector2f, v0::SVector2f, x1::SVector2f, v1::SVector2f)
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
    return t_star, cost(t_star)
end

# see ICRA paper: D.J.Webb et al Kinodynamic RRT* (2013)
"""
function forward_reachable_box(x0::SVector2f, v0::SVector2f, r::Float64)
    @fastmath tau_x_plus = 2/3*(-(v0.^2).+r + v0.*sqrt.((v0.^2).+r))
    @fastmath tau_x_minus = 2/3*(-(v0.^2).+r - v0.*sqrt.((v0.^2).+r))
    @fastmath xmax = v0.*tau_x_plus + x0 + sqrt.(1/3*(tau_x_plus.^2).*(-tau_x_plus.+r))
    @fastmath xmin = v0.*tau_x_minus + x0 - sqrt.(1/3*(tau_x_minus.^2).*(-tau_x_minus.+r))

    @fastmath tau_v_plus = 0.5*r
    @fastmath vmax = v0 .+ sqrt.(tau_v_plus.*(-tau_v_plus.+r))
    @fastmath vmin = v0 .- sqrt.(tau_v_plus.*(-tau_v_plus.+r))
    return xmin, xmax, vmin, vmax
end
"""

function forward_reachable_box(x0::SVector2f, v0::SVector2f, r::Float64)
    @fastmath tau_x_plus = 2/3*(-(v0.^2).+r + v0.*sqrt.((v0.^2).+r))
    @fastmath tau_x_minus = 2/3*(-(v0.^2).+r - v0.*sqrt.((v0.^2).+r))
    @fastmath xmax = v0.*tau_x_plus + x0 + sqrt.(1/3*(tau_x_plus.^2).*(-tau_x_plus.+r))
    @fastmath xmin = v0.*tau_x_minus + x0 - sqrt.(1/3*(tau_x_minus.^2).*(-tau_x_minus.+r))

    @fastmath tau_v_plus = 0.5*r
    @fastmath vmax = v0 .+ sqrt.(tau_v_plus.*(-tau_v_plus.+r))
    @fastmath vmin = v0 .- sqrt.(tau_v_plus.*(-tau_v_plus.+r))
    return xmin, xmax, vmin, vmax
end

function backward_reachable_box(x0::SVector2f, v0::SVector2f, r::Float64)
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
    s_set_filtered = SVector4f[]
    for s in s_set
        @views xq = SVector2f(s[1:2])
        @views vq = SVector2f(s[3:4])
        @views xc = SVector2f(s_c[1:2])
        @views vc = SVector2f(s_c[3:4])
        ans = find_tau_star(xc, vc, xq, vq)
        ans[2]<r && push!(s_set_filtered, s)
    end
    return s_set_filtered
end

function filter_freachable(s_set::Vector{SVector4f}, s_c::SVector4f, r::Float64, ForR::Symbol)
    s_set_filtered = SVector4f[]
    @views x_c = SVector2f(s_c[1:2])
    @views v_c = SVector2f(s_c[3:4])
    xmin, xmax, vmin, vmax = ((ForR==:F) ?
                              forward_reachable_box(x_c, v_c, r)
                              :
                              backward_reachable_box(x_c, v_c, r))
                              
    function isinside(s::SVector4f)
        @inbounds !(xmin[1]<s[1]<xmax[1]) && return false
        @inbounds !(xmin[2]<s[2]<xmax[2]) && return false
        @inbounds !(vmin[1]<s[3]<vmax[1]) && return false
        @inbounds !(vmin[2]<s[4]<vmax[2]) && return false
        return true
    end
    for s in s_set
        isinside(s) && push!(s_set_filtered, s)
    end
    return s_set_filtered
end

function test()
    N = 10^5
    s_set = SVector4f[]
    for i = 1:N
        push!(s_set, SVector4f(rand()-0.5, rand()-0.5, 3*rand()-1.5, 3*rand()-1.5))
    end
    s_c = SVector4f(+0, +0, 0.2, 0.2)
    s_filter = filter_freachable(s_set, s_c, 1.0, :F)
    s_filter2 = filter_freachable(s_set, s_c, 1.0, :B)

    ss = zeros(2, length(s_filter))
    ss2 = zeros(2, length(s_filter2))
    for i = 1:length(s_filter)
        ss[:, i] = s_filter[i][1:2]
    end
    for i = 1:length(s_filter2)
        ss2[:, i] = s_filter2[i][1:2]
    end
    scatter(ss[1, :], ss[2, :], s=3)
    scatter(ss2[1, :], ss2[2, :], s=3)
    xlim(-0.6, 0.6)
    ylim(-0.6, 0.6)
    #return s_filter
end

