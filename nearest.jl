function fuck()
    x_set = rand(2, 10000)
    x0 =  [0, 0]
    N = 10000
    for j = 1:N
        for i=1:N
            xx = x_set[1, i]
            yy = x_set[2, i]
            dist = sqrt((x0[1]-xx)^2+(x0[2]-yy)^2)
        end
        println(j)
    end
end

function hack()
    r = 0.05
    N = 10000
    x_set = Array{Float64, 2}[]
    for i=1:N
        x = rand(2, 1)
        push!(x_set, x)
    end

    set = Int64[]
    for j=1:N
        for i=1:N
            dist = mynorm(x_set[i], x_set[1])
            if dist<r
                push!(set, i)
            end
        end
    end
    
end
@inline mynorm(p, q) = sqrt((p[1]-q[1])^2+(p[2]-q[2])^2)


@time hack()
