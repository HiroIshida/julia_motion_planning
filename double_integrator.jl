const Tuple2f = Tuple{Float64, Float64}
Base.:+(x::Tuple2f, y::Tuple2f) = (x[1]+y[1], x[2]+y[2])
Base.:+(x::Tuple2f, y::Float64) = (x[1]+y, x[2]+y)
Base.:+(y::Float64, x::Tuple2f) = (x[1]+y, x[2]+y)

Base.:-(x::Tuple2f, y::Tuple2f) = (x[1]-y[1], x[2]-y[2])
Base.:-(x::Tuple2f, y::Float64) = (x[1]-y, x[2]-y)
Base.:-(y::Float64, x::Tuple2f) = (y-x[1], y-x[2])

Base.:-(x::Tuple2f) = (-x[1], -x[2])
Base.:/(x::Tuple2f, y::Float64) = (x[1]/y, x[2]/y)
Base.:/(x::Tuple2f, y::Int64) = (x[1]/y, x[2]/y)
Base.:*(x::Tuple2f, y::Float64) = (x[1]*y, x[2]*y)
Base.:*(x::Tuple2f, y::Int64) = (x[1]*y, x[2]*y)
Base.:*(x::Float64, y::Tuple2f) = (x*y[1], x*y[2])
Base.:*(x::Int64, y::Tuple2f) = (x*y[1], x*y[2])
dot(x::Tuple2f, y::Tuple2f) = x[1]*y[1]+x[2]*y[2]

broadcast(*, x::Tuple2f, y::Tuple2f) = (x[1]*y[1], x[2]*y[2])
broadcast(sqrt, x::Tuple2f) = (sqrt(x[1]), sqrt(x[2]))



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

function find_tau_star(x0::Tuple2f, v0::Tuple2f, x1::Tuple2f, v1::Tuple2f)
    x01 = x1 - x0
    v01 = v1 - v0
    p = -4*(dot(v0, v0)+dot(v1, v1)+dot(v0, v1))
    q = 24*dot(v0+v1, x01)
    r = -36*dot(x01, x01)
    cost(t) = t + dot(v01, 4.0*v01/t-6(-v0+x01)/t^2)+dot(-6*v01/t^2+12(x01-v0*t)/t^3, -v0*t+x01)
    @fastmath d_cost(t) = t^4+p*t*t+q*t+r # df(t)/dt
    @fastmath dd_cost(t) = 4.0*t*t*t+2*p*t+q # ddf(t)/dt^2

    t_min = 0.0
    t_max = 10.0
    t_star = bisection_newton(d_cost, dd_cost, t_min, t_max)
    return t_star, cost(t_star)
end

# see ICRA paper: D.J.Webb et al Kinodynamic RRT* (2013)
function forward_reachable_box(x0::Tuple2f, v0::Tuple2f, r::Float64)
    r = 2.0
    x0 = (0.0, 0.0)
    v0 = (0.0, 0.0)
    #tau = 2/2*(r-v0
    tau_x_plus = 2/3*(r-v0.^2 + v0.*sqrt.(r+v0.^2))
    tau_x_minus = 2/3*(r-v0.^2 - v0.*sqrt.(r+v0.^2))
    xmax = v0.*tau_x_plus + x0 + sqrt.(1/3*(tau_x_plus.^2).*(r-tau_x_plus))
    xmin = v0.*tau_x_minus + x0 - sqrt.(1/3*(tau_x_plus.^2).*(r-tau_x_minus))

    tau_v_plus = 0.5*r
    vmax_tmp = v0 + sqrt.(tau_v_plus.*(r-tau_v_plus))
    vmin_tmp = v0 - sqrt.(tau_v_plus.*(r-tau_v_plus))
    vmax = (vmax_tmp, vmax_tmp)
    vmin = (vmin_tmp, vmin_tmp)
    return xmin, xmax, vmin, vmax
end

function is_forwardreachable(xq, vq, x_c, v_c, r)
    # xq, vq (querry)
    #x_c, v_c (start point)
    xmin, xmax, vmin, vmax = forward_reachable_box(xc, vc,r)
    !(x_min[1]<xs[1]<x_max[1]) && return false
    !(x_min[2]<xs[2]<x_max[2]) && return false
    !(v_min[1]<vs[1]<v_mav[1]) && return false
    !(v_min[2]<vs[2]<v_mav[2]) && return false
    return true
end

function fliter_forwardreachable(s_set, s_c, r)
    x_set_filtered = Tuple2f[]
    v_set_filtered = Tuple2f[]
    for s in s_set

    end
end


function hage()
    N = 10^6
    s = 0.0
    for i=1:N
        ans = find_tau_star((0.0, 0.0), (0.0, 0.0), (0.2, 0.2), (0.01, 0.05))
    end
    ans = find_tau_star((0.0, 0.0), (0.0, 0.0), (0.2, 0.2), (0.01, 0.05))
    x0 = (0.0, 0.0)
    v0 = (0.0, 0.0)
    forward_reachable_box(x0, v0, 1.0)
    println(ans)
end

@time hage()



