"""
VectorizationTransformations.jl is a package for working with linear transformations of vectorizations.
"""
module VectorizationTransformations

using LinearAlgebra
using SparseArrays

export vech, duplication_matrix, elimination_matrix, commutation_matrix, symmetrizer_matrix


"""
    vech(A::AbstractMatrix{T}) where {T}

Return a vector containing the elements of the lower triangle (including the diagonal) of the symmetric matrix `A`.

# Arguments
- `A::AbstractMatrix{T}`: A symmetric matrix of any type `T`.

# Returns
- `Vector{T}`: A vector of length `n * (n + 1) ÷ 2`, where `n` is the number of rows (or columns) in `A`. This vector contains the elements from the lower triangle of `A`.


# Examples
```julia
A = [1 2 3; 2 4 5; 3 5 6]
v = vech(A) # Returns [1, 2, 3, 4, 5, 6]
```
"""
function vech(A::AbstractMatrix{T}) where {T}
    # Use issymmetric to check if the matrix is symmetric
    !issymmetric(A) && error("The matrix is not symmetric.")

    n = size(A, 1)
    len = n * (n + 1) ÷ 2
    v = Vector{T}(undef, len)

    k = 1
    @inbounds for j in 1:n
        for i in j:n
            v[k] = A[i, j]
            k += 1
        end
    end

    return v
end


"""
    duplication_matrix([T], n::Int)

Create a sparse matrix used to duplicate the elements of a symmetric matrix of size `n x n`.

This function generates a sparse matrix `D` such that `D * vech(A) = vec(A)`. 

# Arguments
- `T`: The type of the elements in the resulting sparse matrix. Optional, defaults to `Int64`.
- `n::Int`: The size of the symmetric matrix (number of rows or columns).

# Returns
- `SparseMatrixCSC{T,Int64}`: A sparse matrix of size `n^2 x n * (n + 1) ÷ 2`.

See also: `vech`, `elimination_matrix`

# Examples
```julia
A = [1 2 3; 2 4 5; 3 5 6]
D = duplication_matrix(3)
D * (1:6) == vec(A)
```
"""
function duplication_matrix(T, n::Int)
    len = n * (n + 1) ÷ 2
    rows = Int[]
    cols = Int[]
	sizehint!(rows,n^2)
	sizehint!(cols,len)
	
    k = 1
    for i in 1:n
        for j in i:n
            push!(rows, (i - 1) * n + j)
            push!(cols, k)

            if i != j
                push!(rows, (j - 1) * n + i)
                push!(cols, k)
            end

            k += 1
        end
    end

    return sparse(rows, cols, one(T), n^2, len)
end
duplication_matrix(n::Int) = duplication_matrix(Int64,n)




"""
    elimination_matrix([T], n::Int)

Create a sparse matrix used to duplicate the elements of a symmetric matrix of size `n x n`.

This function generates a sparse matrix `L` such that `L * vec(A) = vech(A)`. 

# Arguments
- `T`: The type of the elements in the resulting sparse matrix. Optional, defaults to `Int64`.
- `n::Int`: The size of the symmetric matrix (number of rows or columns).

# Returns
- `SparseMatrixCSC{T,Int64}`: A sparse matrix of size `n * (n + 1) ÷ 2 x n^2`.

See also: `vech`, `duplication_matrix`

# Examples
```julia
A = [1 2 3; 2 4 5; 3 5 6]
L = elimination_matrix(3)
L * vec(A) == 1:6
```
"""
function elimination_matrix(T, n::Int)
    len = n * (n + 1) ÷ 2
    rows = Int[]
    cols = Int[]
	sizehint!(rows,len)
	sizehint!(cols,n^2)
	

    for j in 1:n
        for i in j:n
            push!(rows, (j-1)*n + i - ((j-1)*j)÷2)
            push!(cols, (j-1)*n + i)
        end
    end

    return sparse(rows, cols, one(T), len, n^2)
end
elimination_matrix(n::Int) = elimination_matrix(Int64,n)




"""
    commutation_matrix([T], m::Int, n::Int = m)

Create a sparse commutation matrix `K` such that `K * vec(A) = vec(A')`.

# Arguments
- `T`: The type of the elements in the resulting sparse matrix. Optional, defaults to `Int64`.
- `m::Int`: The number of rows in the matrices to be commuted.
- `n::Int`: The number of columns in the matrices to be commuted. Defaults to `m` if not provided.

# Returns
- `SparseMatrixCSC{T,Int64}`: A sparse matrix of size `m*n x m*n`.

See also: `symmetrizer_matrix`

# Examples
```julia
A = reshape(1:6,3,2)
K = commutation_matrix(3,2)
K * vec(A) == vec(A')
```
"""
function commutation_matrix(T, m::Int, n::Int=m)
    rows = 1:(m*n)
    cols = Int[]
	sizehint!(cols,m*n)
	
    for k in 1:m
        for l in k:m:(m*n)
            push!(cols, l)
        end
    end

    return sparse(rows, cols, one(T), m*n, m*n)
end
commutation_matrix(m::Int, n::Int=m) = commutation_matrix(Int64,m,n)





"""
    symmetrizer_matrix(n::Int)

Create a sparse symmetrizer matrix `S` such that `S * vec(A) = vec((A+A')/2)`.

# Arguments
- `n::Int`: The size of the symmetric matrix (number of rows or columns).


# Returns
- `SparseMatrixCSC{Float64,Int64}`: A sparse matrix of size `n^2 x n^2`.

See also: `commutation_matrix`

# Examples
```julia
A = rand(3,3)
S = symmetrizer_matrix(3)
S * vec(A) == vec((A+A')/2)
```
"""
symmetrizer_matrix(n::Int) = 1/2  * (commutation_matrix(n) + I)




end # module
