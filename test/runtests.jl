using PicardLefschetz
using Test

@testset "PicardLefschetz.jl" begin
    # Write your tests here.

    @test square(3) == 9
    @test square(π) ≈ 9.8696 rtol=1e-4

end
