function Tensor(N::Int, T::DataType)
    N==0 && return eval(T)
    return Tensor(N-1, eval(Expr(:curly, :Vector, T)))
end

data = Tensor(3, Float64)[]
