# Define simple generic types and methods for testing

sprintcompact(x) = sprint(show, x; context=:compact=>true)

struct TestSharpness <: TreatmentSharpness end
struct TestParaCondition <: ParallelCondition end
struct TestParaStrength <: ParallelStrength end

struct TestTreatment <: AbstractTreatment
    time::Symbol
    ref::Int
end

ttreat(time::Term, ref::ConstantTerm) = TestTreatment(time.sym, ref.n)

struct TestParallel{C,S} <: AbstractParallel{C,S}
    e::Int
end

TestParallel(e::Int) = TestParallel{ParallelCondition,ParallelStrength}(e)
tpara(c::ConstantTerm) = TestParallel{ParallelCondition,ParallelStrength}(c.n)

teststep(tr::AbstractTreatment, pr::AbstractParallel) =
    (str=sprint(show, tr), spr=sprint(show, pr))
const TestStep = StatsStep{:TestStep, typeof(teststep), true}
required(::TestStep) = (:tr, :pr)

testnextstep(::AbstractTreatment, str::String) = (next="next"*str,)
const TestNextStep = StatsStep{:TestNextStep, typeof(testnextstep), true}
required(::TestNextStep) = (:tr, :str)

const TestDID = DiffinDiffsEstimator{:TestDID, Tuple{TestStep,TestNextStep}}
const NotImplemented = DiffinDiffsEstimator{:NotImplemented, Tuple{}}

const TR = TestTreatment(:t, 0)
const PR = TestParallel(0)

struct TestResult{TR} <: DIDResult{TR}
    coef::Vector{Float64}
    vcov::Matrix{Float64}
    vce::Nothing
    tr::AbstractTreatment
    nobs::Int
    dof_residual::Int
    yname::String
    coefnames::Vector{String}
    coefinds::Dict{String, Int}
    treatcells::VecColumnTable
    weightname::Symbol
    extra1::Int
    extra2::Symbol
    extra3::String
    extra4::Vector{Float64}
    extra5::Vector{String}
    extra6::Matrix{Float64}
    extra7::Vector{Symbol}
    extra8::Nothing
    extra9::Array{Float64,3}
end

@fieldequal TestResult

function TestResult(n1::Int, n2::Int)
    N = n1*(n2+1)
    coef = collect(Float64, 1:N)
    tr = TR
    tnames = ["rel: $a & c: $b" for a in 1:n1 for b in 1:n2]
    cnames = vcat(tnames, ["c"*string(i) for i in n1*n2+1:N])
    cinds = Dict(cnames .=> 1:N)
    tcells = VecColumnTable((rel=repeat(1:n1, inner=n2), c=repeat(1:n2, outer=n1)))
    return TestResult{typeof(tr)}(coef, coef.*coef', nothing, tr, N, N-1, "y", cnames, cinds,
        tcells, :w, 1, :a, "a", [1.0], ["a"], [1.0 2.0], [:a], nothing, ones(1,1,1))
end

struct TestResultBARE{TR} <: DIDResult{TR}
    coef::Vector{Float64}
    vcov::Matrix{Float64}
    vce::Nothing
    tr::AbstractTreatment
    nobs::Int
    yname::String
    coefnames::Vector{String}
    treatcells::VecColumnTable
    weightname::Symbol
end

@fieldequal TestResultBARE

function TestResultBARE(n1::Int, n2::Int)
    N = n1*n2
    coef = collect(Float64, 1:N)
    tr = dynamic(:time, (-1,))
    tnames = ["rel: $a & c: $b" for a in 1:n1 for b in 1:n2]
    tcells = VecColumnTable((rel=repeat(1:n1, inner=n2), c=repeat(1:n2, outer=n1)))
    return TestResultBARE{typeof(tr)}(coef, coef.*coef', nothing, tr, N, "y", tnames, tcells, :w)
end

function result(::Type{TestDID}, nt::NamedTuple)
    return merge(nt, (result=TestResult(2, 2),))
end
