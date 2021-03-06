@testset "cb" begin
    @test cb() == []

    f = @formula(y ~ cb(1))
    @test f.rhs.forig(f.rhs.args_parsed...) == [1]

    f = @formula(y ~ cb(1,2))
    @test f.rhs.forig(f.rhs.args_parsed...) == [1,2]

    f = @formula(y ~ cb(1,x))
    @test_throws ArgumentError f.rhs.forig(f.rhs.args_parsed...)
end

@testset "unpack" begin
    @test unpack(term(1)) == 1
    @test unpack(term(:x)) == :x
    @test unpack(term(:Unconditional)) == Unconditional()
    f = @formula(y ~ cb)
    @test unpack(f.rhs) == []
    f = @formula(y ~ cb(1,2))
    @test unpack(f.rhs) == [1,2]
end

@testset "exampledata" begin
    @test exampledata() == (:hrs, :nsw, :mpdta)
    hrs = exampledata("hrs")
    @test size(hrs) == (3280, 11)
    @test hrs[1] isa Vector
    @test size(exampledata(:nsw)) == (32834, 11)
    @test size(exampledata(:mpdta)) == (2500, 5)
end

@testset "checktable" begin
    @test_throws ArgumentError checktable(rand(3,3))
    @test_throws ArgumentError checktable([(a=1, b=2), (a=1, b=2)])
end
