bool = trues(10000)
idxset_unvisited = findall(bool)
p = [0 0]
println("fuck")
@time for idx in idxset_unvisited
    dist2(p, this.Pset[idx])
end
println("fuck")
@time for i = [x for x in 1:10000]
    dist2(p, this.Pset[i])
end
