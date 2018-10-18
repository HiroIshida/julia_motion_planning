#import StaticArrays: SVector
using StaticArrays

const vec2d = SVector{2, Float64}
const smat{N, M} = SMatrix{N, M} where {N, M}

struct data
    b_min::vec2d
    b_max::vec2d
    Nx::Int64
    Ny::Int64
    TPM::SMatrix

    function data(b_min::vec2d, b_max::vec2d, Nx::Int64, Ny::Int64)
        TPM = compute_TPM(b_min, b_max, Nx, Ny)
        #TPM = @SMatrix randn(Nx, Ny);
        new(b_min, b_max, Nx, Ny, TPM)
    end
end

function compute_TPM(b_min, b_max, Nx, Ny)

end



b_min = vec2d(0, 0)
b_max = vec2d(10, 10)
Nx = 10
Ny = 10
D=data(b_min, b_max, Nx, Ny)





