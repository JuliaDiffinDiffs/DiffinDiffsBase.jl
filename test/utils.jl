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
    @test size(exampledata(:hrs)) == (3280,)
    @test size(exampledata(:nsw)) == (32834,)
    @test size(exampledata(:mpdta)) == (2500,)
end

@testset "RotatingTimeValue" begin
    @test rotatingtime(1.0, Date(1)) == RotatingTimeValue(1.0, Date(1))
    rot = repeat(5:-1:1, inner=5)
    time = repeat(1:5, outer=5)
    rt = rotatingtime(rot, time)
    @test eltype(rt) == RotatingTimeValue{Int, Int}
    @test getfield.(rt, :rotation) == rot
    @test getfield.(rt, :time) == time

    @test rt[1] - rt[2] == -1
    @test_throws ArgumentError rt[1] - rt[6]

    @test isless(rt[1], rt[2])
    @test isless(rt[6], rt[1])
    @test rt[1] == RotatingTimeValue(Int32(5), Int32(1))

    @test sprint(show, rt[1]) == "5_1"
    w = VERSION < v"1.6.0" ? "" : " "
    @test sprint(show, MIME("text/plain"), rt[1]) == """
        RotatingTimeValue{Int64,$(w)Int64}:
          rotation: 5
          time:     1"""
end