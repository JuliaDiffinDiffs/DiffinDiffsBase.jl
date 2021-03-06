@testset "_mult!" begin
    a = collect(0:5)
    b = collect(1:6)
    _mult!(a, b, 6)
    @test a == 0:7:35
end

@testset "findcell" begin
    hrs = exampledata("hrs")
    @test_throws ArgumentError findcell((), hrs)
    @test_throws ArgumentError findcell((:wave,), hrs, falses(size(hrs, 1)))

    rows = findcell((:wave,), hrs)
    @test length(rows) == 5
    @test rows[UInt32(5)] == findall(x->x==7, hrs.wave)

    df = DataFrame(hrs)
    df.wave = PooledArray(df.wave)
    @test findcell((:wave,), df) == rows

    esample = hrs.wave.!=11
    rows = findcell((:wave, :wave_hosp), hrs, esample)
    @test length(rows) == 16
    @test rows[one(UInt32)] == intersect(findall(x->x==10, view(hrs.wave, esample)), findall(x->x==10, view(hrs.wave_hosp, esample)))

    rows = findcell((:wave, :wave_hosp, :male), hrs)
    @test length(rows) == 40
end

@testset "cellrows" begin
    hrs = exampledata("hrs")
    cols0 = subcolumns(hrs, (:wave, :wave_hosp), falses(size(hrs, 1)))
    cols = subcolumns(hrs, (:wave, :wave_hosp))
    rows_dict0 = IdDict{UInt32, Vector{Int}}()
    @test_throws ArgumentError cellrows(cols, rows_dict0)
    rows_dict = findcell(cols)
    @test_throws ArgumentError cellrows(cols0, rows_dict)

    cells, rows = cellrows(cols, rows_dict)
    @test length(cells[1]) == length(rows) == 20
    @test Tables.matrix(cells) ==
        sortslices(unique(hcat(hrs.wave, hrs.wave_hosp), dims=1), dims=1)
    @test propertynames(cells) == [:wave, :wave_hosp]
    @test rows[1] == intersect(findall(x->x==7, hrs.wave), findall(x->x==8, hrs.wave_hosp))

    df = DataFrame(hrs)
    df.wave = ScaledArray(hrs.wave, 1)
    df.wave_hosp = ScaledArray(hrs.wave_hosp, 1)
    cols1 = subcolumns(df, (:wave, :wave_hosp))
    cells1, rows1 = cellrows(cols1, rows_dict)
    @test rows1 == rows
    @test cells1 == cells
    @test cells1.wave isa ScaledArray
    @test cells1.wave_hosp isa ScaledArray

    rot = ones(size(df, 1))
    df.wave = RotatingTimeArray(rot, hrs.wave)
    df.wave_hosp = RotatingTimeArray(rot, hrs.wave_hosp)
    cols2 = subcolumns(df, (:wave, :wave_hosp))
    cells2, rows2 = cellrows(cols2, rows_dict)
    @test rows2 == rows
    @test cells2 == cells
    @test cells2.wave isa RotatingTimeArray
    @test cells2.wave_hosp isa RotatingTimeArray
end

@testset "settime aligntime" begin
    hrs = exampledata("hrs")
    t = settime(hrs.wave, 2, reftype=Int16)
    @test eltype(refarray(t)) == Int16
    @test sort!(unique(refarray(t))) == 1:3
    t = settime(hrs.wave, 0.5)
    @test eltype(t) == Float64

    rot = isodd.(hrs.wave)
    t = settime(hrs.wave, rotation=rot)
    @test eltype(t) == RotatingTimeValue{Bool, eltype(hrs.wave)}
    @test refpool(t) === nothing
    @test eltype(refarray(t.time)) == Int32
    @test t.rotation == rot
    @test t.time == hrs.wave

    df = DataFrame(hrs)
    df.wave1 = Date.(df.wave)
    t = settime(df.wave1, Year(1))
    @test t == df.wave1
    @test t.pool == Date(7):Year(1):Date(11)

    t = settime(df.wave1, Year(1), rotation=rot)
    @test eltype(t) == RotatingTimeValue{Bool, eltype(df.wave1)}
    @test eltype(refarray(t.time)) == Int32
    @test t.time isa ScaledArray

    df.t1 = settime(df.wave, 1)
    df.t2 = settime(df.wave, 2)
    t = aligntime(df, :t2, :t1)
    @test t == df.t2
    @test t.pool == df.t1.pool
    @test t.invpool == df.t1.invpool
    df.t2 = settime(df.wave, start=0, stop=20)
    t = aligntime(df, :t2, :t1)
    @test t == df.t2
    @test t.pool == df.t1.pool

    t1 = settime(df.wave1, Year(1), rotation=rot)
    t2 = settime(df.wave1, Year(1), start=Date(5))
    t = aligntime(t2, t1)
    @test t == t1
    @test t.time.pool == t1.time.pool
    t2 = settime(df.wave1, Year(1), start=Date(5), rotation=rot)
    t = aligntime(t2, t1)
    @test t == t1
    @test t.time.pool == t1.time.pool
end

@testset "PanelStructure" begin
    hrs = exampledata("hrs")
    N = size(hrs, 1)
    panel = setpanel(hrs, :hhidpn, :wave)
    @test length(unique(panel.refs)) == N
    inidbounds = view(hrs.hhidpn, 2:N) .== view(hrs.hhidpn, 1:N-1)
    @test view(diff(panel.refs), inidbounds) == view(diff(hrs.wave), inidbounds)
    @test length(panel.idpool) == 656
    @test panel.timepool == 7:11

    panel1 = setpanel(hrs, :hhidpn, :wave, step=0.5, reftype=Int)
    @test view(diff(panel1.refs), inidbounds) == 2 .* view(diff(hrs.wave), inidbounds)
    @test panel1.timepool == 7.0:0.5:11.0
    @test eltype(panel1.refs) == Int

    @test_throws ArgumentError setpanel(hrs, :hhidpn, :oop_spend)
    @test_throws DimensionMismatch setpanel(hrs.hhidpn, 1:100)

    @test sprint(show, panel) == "Panel Structure"
    t = VERSION < v"1.6.0" ? "Array{Int64,1}" : " Vector{Int64}"
    @test sprint(show, MIME("text/plain"), panel) == """
        Panel Structure:
          idpool:   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  647, 648, 649, 650, 651, 652, 653, 654, 655, 656]
          timepool:   7:1:11
          laginds:  Dict{Int64,$t}()"""

    lags = findlag!(panel)
    @test sprint(show, MIME("text/plain"), panel) == """
        Panel Structure:
          idpool:   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  647, 648, 649, 650, 651, 652, 653, 654, 655, 656]
          timepool:   7:1:11
          laginds:  Dict(1 => [4, 5, 1, 2, 0, 7, 0, 10, 8, 6  …  0, 3275, 3271, 3273, 3274, 3279, 3280, 3277, 3278, 0])"""

    leads = findlead!(panel, -1)
    @test leads == leads
    @test_throws ArgumentError findlag!(panel, 5)

    @test ilag!(panel) === panel.laginds[1]
    @test ilag!(panel, 2) === panel.laginds[2]
    @test ilead!(panel) === panel.laginds[-1]

    l1 = lag(panel, hrs.wave)
    @test view(l1, hrs.wave.!=7) == view(hrs.wave, hrs.wave.!=7) .- 1
    @test all(ismissing, view(l1, hrs.wave.==7))
    l2 = lag(panel, hrs.wave, 2, default=-1)
    @test view(l2, hrs.wave.>=9) == view(hrs.wave, hrs.wave.>=9) .- 2
    @test all(x->x==-1, view(l2, hrs.wave.<9))
    @test eltype(l2) == Int
    l_2 = lead(panel, hrs.wave, 2, default=-1.0)
    @test view(l_2, hrs.wave.<=9) == view(hrs.wave, hrs.wave.<=9) .+ 2
    @test all(x->x==-1, view(l_2, hrs.wave.>9))
    @test eltype(l_2) == Int

    @test_throws DimensionMismatch lag(panel, 1:10)

    d1 = diff(panel, hrs.wave)
    @test all(x->x==1, view(d1, hrs.wave.!=7))
    @test all(ismissing, view(d1, hrs.wave.==7))
    d2 = diff(panel, hrs.wave, l=2, default=-1)
    @test all(x->x==2, view(d2, hrs.wave.>=9))
    @test all(x->x==-1, view(d2, hrs.wave.<9))
    @test eltype(d2) == Int
    d_2 = diff(panel, hrs.wave, l=-2, default=-1.0)
    @test all(x->x==-2, view(d_2, hrs.wave.<=9))
    @test all(x->x==-1, view(d_2, hrs.wave.>9))
    @test eltype(d_2) == Int
    
    d2 = diff(panel, hrs.wave, order=2)
    @test all(x->x==0, view(d2, hrs.wave.>=9))
    @test all(ismissing, view(d2, hrs.wave.<9))
    d2 = diff(panel, hrs.wave, order=2, l=2)
    @test all(x->x==0, view(d2, hrs.wave.==11))
    @test all(ismissing, view(d2, hrs.wave.!=11))

    @test_throws ArgumentError diff(panel, hrs.wave, l=5)
    @test_throws ArgumentError diff(panel, hrs.wave, order=5)
end
