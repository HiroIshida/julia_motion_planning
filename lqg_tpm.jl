mutable struct LQG_TPM
    lo
    hi
    N_split
    dim_state
    TPM
    function LQG_TPM(A, B, Q, R
                     lo = [-1, -1],
                     hi = [1, 1],
                     N_split = [20, 20],
                     dim_state = 2)
        new(lo, hi, N_split, dim_state)
    end
    
end
function y = 


end
