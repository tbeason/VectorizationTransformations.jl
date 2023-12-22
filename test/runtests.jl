using VectorizationTransformations
using Test
using LinearAlgebra

# write tests here

sym(A) = (A + A') / 2

n = 5
A = rand(n,n)


@test duplication_matrix(n)*vech(sym(A)) == vec(sym(A))

@test elimination_matrix(n)*vec(sym(A)) == vech(sym(A))

@test elimination_matrix(n)*duplication_matrix(n) == I(n * (n + 1) ÷ 2)

@test commutation_matrix(n)*vec(A) == vec(A')

@test symmetrizer_matrix(n)*vec(A) == vec(sym(A))

@test reshape(symmetrizer_matrix(n)*vec(A),n,n) == sym(A)

@test  elimination_matrix(n) * symmetrizer_matrix(n) * vec(A) == vech(sym(A))



B = reshape(1:6,3,2)

@test_throws "The matrix is not symmetric." vech(B)

@test commutation_matrix(3,2)*vec(B) == vec(B')

@test commutation_matrix(2,3)*vec(B') == 1:6

@test commutation_matrix(3,2)' * commutation_matrix(3,2) == I

@test commutation_matrix(3,2) * commutation_matrix(3,2)' == I





























## NOTE add JET to the test environment, then uncomment
# using JET
# @testset "static analysis with JET.jl" begin
#     @test isempty(JET.get_reports(report_package(VectorizationTransformations, target_modules=(VectorizationTransformations,))))
# end

## NOTE add Aqua to the test environment, then uncomment
# @testset "QA with Aqua" begin
#     import Aqua
#     Aqua.test_all(VectorizationTransformations; ambiguities = false)
#     # testing separately, cf https://github.com/JuliaTesting/Aqua.jl/issues/77
#     Aqua.test_ambiguities(VectorizationTransformations)
# end
