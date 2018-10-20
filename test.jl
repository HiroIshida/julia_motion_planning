function test()
    p = 1.0
    q = 2.0
    r = 3.0

    function df(t)
        return t^4+p*t*t+q*t+r # df(t)/dt
    end

    function ddf(t)
        return 4.0*t^3+2*p*t+q # ddf(t)/dt^2
    end

    s = 0.0
    t = 1.0
    for i=1:10*6
        s+=df(t)/ddf(t)
        s+=df(t)
    end
    return s
end

@time test()
