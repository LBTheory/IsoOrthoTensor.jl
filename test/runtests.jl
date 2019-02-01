#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

using Test
using IsoOrthoTensor


#------------------------------------------------------------------------------#
#                            Isotropic Tensor Tests                            #
#------------------------------------------------------------------------------#

@testset "Isotropic Tensor Tests                            " begin
    @test Δ(1) == K(2)
    @test Δ(2) == 𝕔((K(2), K(2)), (1,))
    @test Δ(3) == 𝕔((K(2), Δ(2)), (1,))
    @test Δ(4) == 𝕔((K(2), Δ(3)), (1,))
    @test Δ(5) == 𝕔((K(2), Δ(4)), (1,))
end


#------------------------------------------------------------------------------#
#                          Orthogonality Tensor Tests                          #
#------------------------------------------------------------------------------#

@testset "Orthogonality Tensor Tests                        " begin
    @test Ο(1) == K(2)
    @test Ο(2) == 𝕡((K(2), K(2)), (1, 3))
    @test Ο(3) == 𝕡((K(2), K(2), K(2)), (1, 3, 5))
    @test Ο(4) == 𝕡((K(2), K(2), K(2), K(2)), (1, 3, 5, 7))
    @test Ο(5) == 𝕡((K(2), K(2), K(2), K(2), K(2)), (1, 3, 5, 7, 9))
end


