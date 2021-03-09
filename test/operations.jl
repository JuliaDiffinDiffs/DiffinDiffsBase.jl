@testset "findcell" begin
    hrs = exampledata("hrs")
    @test_throws ArgumentError findcell((), hrs)
    @test findcell((:wave,), hrs, falses(size(hrs, 1))) ==
        (table(Matrix{Any}(undef, 0, 1), header=[:wave]), Vector{Int}[Int[]])

    cells, rows = findcell((:wave,), hrs)
    @test getfield(cells, :matrix) == reshape(collect(Any, 7:11), 5, 1)
    @test cells.wave == collect(Any, 7:11)
    @test length(rows) == 5
    @test rows[1] == findall(x->x==7, hrs.wave)

    df = DataFrame(hrs)
    df.wave = PooledArray(df.wave)
    @test findcell((:wave,), df) == (cells, rows)

    esample = hrs.wave.!=11
    cells, rows = findcell((:wave, :wave_hosp), hrs, esample)
    @test length(cells) == length(rows) == 16
    @test getfield(cells, :matrix)[1,:] == [7, 8]
    @test rows[1] == intersect(findall(x->x==7, view(hrs.wave, esample)), findall(x->x==8, view(hrs.wave_hosp, esample)))

    df.wave_hosp = PooledArray(df.wave_hosp)
    @test findcell((:wave, :wave_hosp), df, esample) == (cells, rows)

    cells, rows = findcell((:wave, :wave_hosp, :male), hrs)
    @test length(cells) == 40
end
