const Tuple2f = Tuple{Float64, Float64}
Base.:+(x::Tuple2f, y::Tuple2f) = (x[1]+y[1], x[2]+y[2])
Base.:-(x::Tuple2f, y::Tuple2f) = (x[1]-y[1], x[2]-y[2])
Base.:-(x::Tuple2f) = (-x[1], -x[2])
Base.:/(x::Tuple2f, y::Float64) = (x[1]/y, x[2]/y)
Base.:/(x::Tuple2f, y::Int64) = (x[1]/y, x[2]/y)
Base.:*(x::Tuple2f, y::Float64) = (x[1]*y, x[2]*y)
Base.:*(x::Tuple2f, y::Int64) = (x[1]*y, x[2]*y)
Base.:*(x::Float64, y::Tuple2f) = (x*y[1], x*y[2])
Base.:*(x::Int64, y::Tuple2f) = (x*y[1], x*y[2])
dot(x::Tuple2f, y::Tuple2f) = x[1]*y[1]+x[2]*y[2]

function find_tau_star(x0::Tuple2f, v0::Tuple2f, x1::Tuple2f, v1::Tuple2f)
    x01 = x1 - x0
    v01 = v1 - v0
    p = -4*(dot(v0, v0)+dot(v1, v1)+dot(v0, v1))
    q = 24*dot(v0+v1, x01)
    r = -36*dot(x01, x01)

    # f(t) correspoding to cost(t)
    f(t) = t + dot(v01, 4*v01/t-6(-v0+x01)/t^2)+dot(-6*v01/t^2+12(x01-v0*t)/t^3, -v0*t+x01)
    @fastmath df(t) = t^4+p*t*t+q*t+r # df(t)/dt
    @fastmath ddf(t) = 4.0*t*t*t+2*p*t+q # ddf(t)/dt^2

    ## lets find the root of df(t)/dt, which minimize f(t)
    # solve it by combination of bisection and Newton
    eps = 0.05
    right = 5
    left = 0.01
    est = 0.5*(right + left) # initial estimation for the root
    for itr = 1:20
        df_ddf = df(est)/ddf(est)
        if (right-(est-df_ddf))*((est-df_ddf)-left)<0.0 || abs(df_ddf)<(right-left)/4.0
            if df(est)>0
                right = est
            else
                left = est
            end
            df_ddf = est - (right + left)*0.5
        end
        est -= df_ddf
        abs(df_ddf)<eps && return est, f(est)
    end
    return est, f(est)
end

function func()
    N = 1000000
    s = 0.0
    for i=1:N
        ans = find_tau_star((0.0, 0.0), (0.5, 0.5), (0.2, 0.2), (-0.2, 0.0))
        s += ans[1]
    end
    return s
end

@time func()


