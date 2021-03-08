@testset "_relabel!" begin
    hrs = exampledata("hrs")
    refs, invpool, pool = _label(hrs.wave)
    mult = _relabel!(refs, length(pool), refs)
    @test length(unique(refs)) == 5
    @test mult == 25

    mult = 5
    mult = _relabel!(refs, mult, hrs.wave_hosp)
    @test length(unique(refs)) == 20
    @test mult == 20

    @test_throws DimensionMismatch _relabel!(refs, 20, reshape(hrs.wave_hosp, size(hrs, 1), 1))
end

@testset "findcell" begin
    hrs = exampledata("hrs")
    df = DataFrame(hrs)
    nrow = size(df, 1)
    @test_throws ArgumentError findcell((:wave,), hrs, falses(nrow))

    cellrows = findcell((:wave,), hrs)
    cells = sort(collect(keys(cellrows)))
    @test cells == Tuple[(7,), (8,), (9,), (10,), (11,)]
    @test cellrows[(7,)] == findall(x->x==7, hrs.wave)

    df.wave = PooledArray(df.wave)
    @test findcell((:wave,), df) == cellrows

    esample = hrs.wave.!=11
    cellrows = findcell((:wave, :wave_hosp), hrs, esample)
    cells = sort(collect(keys(cellrows)))
    @test length(cells) == 16
    @test cells[1] == (7,8)
    @test cellrows[(7,8)] == intersect(findall(x->x==7, view(hrs.wave, esample)), findall(x->x==8, view(hrs.wave_hosp, esample)))

    df.wave_hosp = PooledArray(df.wave_hosp)
    @test findcell((:wave, :wave_hosp), df, esample) == cellrows
end
