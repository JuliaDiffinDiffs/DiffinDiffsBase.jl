module DiffinDiffsBase

using Base: @propagate_inbounds
using CSV
using CodecZlib: GzipDecompressorStream
using Combinatorics: combinations
using DataAPI
using DataAPI: refarray, refpool, invrefpool
using Dates: Period, TimeType
using LinearAlgebra: Diagonal
using MacroTools: @capture, isexpr, postwalk
using Missings: allowmissing, disallowmissing
using PooledArrays: _label
using Reexport
using StatsBase: CoefTable, Weights, stderror, uweights
using StatsFuns: tdistccdf, tdistinvcdf
@reexport using StatsModels
using StatsModels: Schema
using Tables
using Tables: AbstractColumns, istable, columnnames, getcolumn

import Base: ==, +, -, *, isless, show, parent, view, diff
import Base: eltype, firstindex, lastindex, getindex, iterate, length, sym_in
import StatsBase: coef, vcov, confint, nobs, dof_residual, responsename, coefnames, weights,
    coeftable
import StatsModels: concrete_term, schema, termvars, lag, lead

# Reexport objects from StatsBase
export coef, vcov, stderror, confint, nobs, dof_residual, responsename, coefnames, weights,
    coeftable

export cb,
       ≊,
       exampledata,
       RotatingTimeValue,
       rotatingtime,

       VecColumnTable,
       VecColsRow,
       subcolumns,
       apply,
       apply_and!,
       apply_and,
       TableIndexedMatrix,

       ScaledArray,
       ScaledVector,
       ScaledMatrix,

       TreatmentSharpness,
       SharpDesign,
       sharp,
       AbstractTreatment,
       DynamicTreatment,
       dynamic,

       ParallelCondition,
       Unconditional,
       unconditional,
       CovariateConditional,
       ParallelStrength,
       Exact,
       exact,
       Approximate,
       AbstractParallel,
       TrendParallel,
       NeverTreatedParallel,
       nevertreated,
       NotYetTreatedParallel,
       notyettreated,
       istreated,

       TermSet,
       termset,
       eachterm,
       TreatmentTerm,
       treat,

       findcell,
       cellrows,
       settime,
       PanelStructure,
       setpanel,
       findlag!,
       findlead!,
       ilag!,
       ilead!,
       diff!,
       diff,

       StatsStep,
       AbstractStatsProcedure,
       SharedStatsStep,
       PooledStatsProcedure,
       StatsSpec,
       proceed,
       @specset,

       CheckData,
       GroupTerms,
       CheckVars,
       MakeWeights,

       DiffinDiffsEstimator,
       DefaultDID,
       did,
       didspec,
       @did,
       AbstractDIDResult,
       DIDResult,
       AggregatedDIDResult,
       vce,
       outcomename,
       treatnames,
       treatcells,
       ntreatcoef,
       treatcoef,
       treatvcov,
       coefinds,
       ncovariate,
       agg,
       SubDIDResult,
       TransformedDIDResult,
       TransSubDIDResult,
       lincom,
       rescale

include("utils.jl")
include("tables.jl")
include("ScaledArrays.jl")
include("treatments.jl")
include("parallels.jl")
include("terms.jl")
include("operations.jl")
include("StatsProcedures.jl")
include("procedures.jl")
include("did.jl")

end
