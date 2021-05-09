"""
    checkdata!(args...)

Check `data` is `Tables.AbstractColumns`-compatible
and find valid rows for options `subset` and `weightname`.
See also [`CheckData`](@ref).
"""
function checkdata!(data, subset::Union{BitVector, Nothing}, weightname::Union{Symbol, Nothing})
    checktable(data)
    nrow = Tables.rowcount(data)
    if subset !== nothing
        length(subset) == nrow || throw(DimensionMismatch(
            "data contain $(nrow) rows while subset has $(length(subset)) elements"))
        esample = subset
    else
        esample = trues(nrow)
    end

    # A cache that makes updating BitVector (esample or tr_rows) faster
    # See https://github.com/JuliaData/DataFrames.jl/pull/2726
    aux = BitVector(undef, nrow)

    if weightname !== nothing
        colweights = getcolumn(data, weightname)
        aux .= colweights .> 0
        esample .&= aux
        if Missing <: eltype(colweights)
            aux .= .!ismissing.(colweights)
            esample .&= aux
        end
    end
    sum(esample) == 0 && error("no nonmissing data")
    return (esample=esample, aux=aux)
end

"""
    CheckData <: StatsStep

Call [`DiffinDiffsBase.checkdata!`](@ref)
for some preliminary checks of the input data.
"""
const CheckData = StatsStep{:CheckData, typeof(checkdata!), true}

required(::CheckData) = (:data,)
default(::CheckData) = (subset=nothing, weightname=nothing)

"""
    groupterms(args...)

Return the arguments for allowing later comparisons based on object-id.
See also [`GroupTerms`](@ref).
"""
groupterms(treatintterms::TermSet, xterms::TermSet) =
    (treatintterms = treatintterms, xterms = xterms)

"""
    GroupTerms <: StatsStep

Call [`DiffinDiffsBase.groupterms`](@ref)
to obtain one of the instances of `treatintterms` and `xterms`
that have been grouped by `==`
for allowing later comparisons based on object-id.

This step is only useful when working with [`@specset`](@ref) and [`proceed`](@ref).
"""
const GroupTerms = StatsStep{:GroupTerms, typeof(groupterms), false}

required(::GroupTerms) = (:treatintterms, :xterms)

function checktreatvars(::DynamicTreatment{SharpDesign}, pr::TrendParallel{Unconditional},
        treatvars::Vector{Symbol}, data)
    # treatvars should be cohort and time variables
    col1 = getcolumn(data, treatvars[1])
    col2 = getcolumn(data, treatvars[2])
    T1 = nonmissingtype(eltype(col1))
    T2 = nonmissingtype(eltype(col2))
    T1 == T2 || throw(ArgumentError(
        "nonmissing elements from columns $(treatvars[1]) and $(treatvars[2]) have different types $T1 and $T2"))
    T1 <: Union{Integer, RotatingTimeValue{<:Any, <:Integer}} ||
        col1 isa ScaledArrOrSub && col2 isa ScaledArrOrSub ||
        throw(ArgumentError("columns $(treatvars[1]) and $(treatvars[2]) must either have integer elements or be ScaledArrays; see settime and aligntime"))
    T1 <: ValidTimeType ||
        throw(ArgumentError("column $(treatvars[1]) has unaccepted element type $(T1)"))
    eltype(pr.e) == T1 || throw(ArgumentError("element type $(eltype(pr.e)) of control cohorts from $pr does not match element type $T1 from data; expect $T1"))
    if col1 isa ScaledArrOrSub
        first(DataAPI.refpool(col1)) == first(DataAPI.refpool(col2)) &&
            scale(col1) == scale(col2) || throw(ArgumentError(
                "time values in columns $(treatvars[1]) and $(treatvars[2]) are not aligned; see aligntime"))
    end
end

function _overlaptime(tr::DynamicTreatment, tr_rows::BitVector, data)
    control_time = Set(view(refarray(getcolumn(data, tr.time)), .!tr_rows))
    treated_time = Set(view(refarray(getcolumn(data, tr.time)), tr_rows))
    return intersect(control_time, treated_time), control_time, treated_time
end

function overlap!(esample::BitVector, tr_rows::BitVector, aux::BitVector, tr::DynamicTreatment,
        ::NeverTreatedParallel{Unconditional}, treatname::Symbol, data)
    overlap_time, control_time, treated_time = _overlaptime(tr, tr_rows, data)
    if !(length(control_time)==length(treated_time)==length(overlap_time))
        aux[esample] .= view(refarray(getcolumn(data, tr.time)), esample) .∈ (overlap_time,)
        esample[esample] .&= view(aux, esample)
    end
    tr_rows .&= esample
end

function overlap!(esample::BitVector, tr_rows::BitVector, aux::BitVector, tr::DynamicTreatment,
        pr::NotYetTreatedParallel{Unconditional}, treatname::Symbol, data)
    overlap_time, _c, _t = _overlaptime(tr, tr_rows, data)
    timetype = eltype(overlap_time)
    invpool = invrefpool(getcolumn(data, tr.time))
    e = invpool === nothing ? Set(pr.e) : Set(invpool[c] for c in pr.e)
    if !(timetype <: RotatingTimeValue)
        ecut = invpool === nothing ? pr.ecut[1] : invpool[pr.ecut[1]]
        filter!(x -> x < ecut, overlap_time)
        isvalidcohort = x -> x < ecut || x in e
    else
        ecut = invpool === nothing ? pr.ecut : (invpool[e] for e in pr.ecut)
        ecut = IdDict(e.rotation=>e.time for e in ecut)
        filter!(x -> x.time < ecut[x.rotation], overlap_time)
        isvalidcohort = x -> x.time < ecut[x.rotation] || x in e
    end
    aux[esample] .= view(refarray(getcolumn(data, tr.time)), esample) .∈ (overlap_time,)
    esample[esample] .&= view(aux, esample)
    aux[esample] .= isvalidcohort.(view(refarray(getcolumn(data, treatname)), esample))
    esample[esample] .&= view(aux, esample)
    tr_rows .&= esample
end

"""
    checkvars!(args...)

Exclude rows with missing data or violate the overlap condition
and find rows with data from treated units.
See also [`CheckVars`](@ref).
"""
function checkvars!(data, tr::AbstractTreatment, pr::AbstractParallel,
        yterm::AbstractTerm, treatname::Symbol, esample::BitVector, aux::BitVector,
        treatintterms::TermSet, xterms::TermSet)
    # Do not check eltype of treatintterms
    treatvars = union([treatname], termvars(tr), termvars(pr))
    checktreatvars(tr, pr, treatvars, data)

    allvars = union(treatvars, termvars(yterm), termvars(xterms))
    for v in allvars
        col = getcolumn(data, v)
        if Missing <: eltype(col)
            aux .= .!ismissing.(col)
            esample .&= aux
        end
    end
    # Values of treatintterms from untreated units are ignored
    tr_rows = copy(esample)
    aux[esample] .= istreated.(Ref(pr), view(getcolumn(data, treatname), esample))
    tr_rows[esample] .&= view(aux, esample)
    treatintvars = termvars(treatintterms)
    for v in treatintvars
        col = getcolumn(data, v)
        if Missing <: eltype(col)
            aux[tr_rows] .= .!ismissing.(view(col, tr_rows))
            esample[tr_rows] .&= view(aux, tr_rows)
        end
    end
    isempty(treatintvars) || (tr_rows[tr_rows] .&= view(esample, tr_rows))

    overlap!(esample, tr_rows, aux, tr, pr, treatname, data)
    sum(esample) == 0 && error("no nonmissing data")
    return (esample=esample, tr_rows=tr_rows::BitVector)
end

"""
    CheckVars <: StatsStep

Call [`DiffinDiffsBase.checkvars!`](@ref) to exclude invalid rows for relevant variables.
"""
const CheckVars = StatsStep{:CheckVars, typeof(checkvars!), true}

required(::CheckVars) = (:data, :tr, :pr, :yterm, :treatname, :esample, :aux)
default(::CheckVars) = (treatintterms=TermSet(), xterms=TermSet())
copyargs(::CheckVars) = (6,)

"""
    makeweights(args...)

Construct a generic `Weights` vector.
See also [`MakeWeights`](@ref).
"""
function makeweights(data, esample::BitVector, weightname::Symbol)
    weights = Weights(convert(Vector{Float64}, view(getcolumn(data, weightname), esample)))
    all(isfinite, weights) || error("data column $weightname contain not-a-number values")
    return (weights=weights,)
end

function makeweights(data, esample::BitVector, weightname::Nothing)
    weights = uweights(sum(esample))
    return (weights=weights,)
end

"""
    MakeWeights <: StatsStep

Call [`DiffinDiffsBase.makeweights`](@ref) to create a generic `Weights` vector.
"""
const MakeWeights = StatsStep{:MakeWeights, typeof(makeweights), true}

required(::MakeWeights) = (:data, :esample)
default(::MakeWeights) = (weightname=nothing,)
