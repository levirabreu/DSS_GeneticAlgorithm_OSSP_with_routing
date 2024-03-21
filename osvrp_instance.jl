# Func to read files
using PyCall, BrkgaMpIpr
np = pyimport("numpy")

struct OSVRPInstance <: AbstractInstance

    m::Int
    n::Int
    ϕ::Int

    ρ::Array{Int, 2}
    δ::Array{Int, 2}

    σ::Array{Int, 1}

    O::Array{Tuple{Int, Int}}

    function OSVRPInstance(p_file::String)

        shape = np.loadtxt(p_file, skiprows=1, max_rows=1, dtype="int64")
        m = shape[1]
        n = shape[2]

        ρ = np.loadtxt(p_file, skiprows=4, max_rows=m, dtype="int64")

        σ = np.loadtxt(p_file, skiprows=4 + n + 2, max_rows=1, dtype="int64")

        aux = np.loadtxt(p_file, skiprows=4 + n + 2 + 3, max_rows=1, dtype="int64")
        ϕ = first(aux)

        δ = np.loadtxt(p_file, skiprows=4 + n + 2 + 3 + 3, max_rows=n+1, dtype="int64")

        O = [(i, j) for i in 1:m for j in 1:n]

        new(m, n, ϕ, ρ, δ, σ, O)
    end

end

#instance = OSVRPInstance("Instâncias\\GP03-01.txt")
