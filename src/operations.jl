# Obtain unique labels for row-wise pairs of values from a1 and a2 when mult1 is large enough
function _mult!(a1::Array, mult1::Integer, a2::AbstractArray)
    a1 .+= mult1 .* (a2 .- 1)
end

# A variant of SplitApplyCombine.groupfind using IdDict instead of Dictionaries.Dictionary
function _groupfind(container)
    T = keytype(container)
    inds = IdDict{eltype(container), Vector{T}}()
    @inbounds for i in keys(container)
        push!(get!(Vector{T}, inds, container[i]), i)
    end
    return inds
end

function _cellrows(cols::VecColumnTable, refrows::IdDict)
    ncol = length(cols)
    ncell = length(refrows)
    rows = Vector{Vector{Int}}(undef, ncell)
    cache = Matrix{Any}(undef, ncell, ncol+1)
    r = 0
    @inbounds for (k, v) in refrows
        r += 1
        cache[r, end] = k
        row1 = v[1]
        for c in 1:ncol
            cache[r, c] = cols[c][row1]
        end
    end
    sorted = sortslices(cache, dims=1)
    cells = table(sorted[:,1:ncol], header=columnnames(cols))
    # Collect rows in the same order as cells
    @inbounds for i in 1:ncell
        rows[i] = refrows[sorted[i,end]]
    end
    return cells, rows
end

"""
    findcell(cols::VecColumnTable)
    findcell(names, data, esample=Colon())

Group the row indices of a collection of data columns
so that the row-wise combinations of values from these columns
are the same within each group.

Instead of directly providing the relevant portions of columns as
[`VecColumnTable`](@ref)``,
one may specify the `names` of columns from
`data` of any Tables.jl-compatible table type
over selected rows indicated by `esample`.
Note that unless `esample` covers all rows of `data`,
the row indices are those for the subsample selected based on `esample`
rather than those for the full `data`.

# Returns
- `cells::MatrixTable`: unique row-wise combinations of values from columns.
- `rows::Vector{Vector{Int}}`: row indices for each combination.
"""
function findcell(cols::VecColumnTable)
    ncol = size(cols, 2)
    ncol == 0 && throw(ArgumentError("no data column is found"))
    if size(cols, 1) == 0
        cells = table(Matrix{Any}(undef, 0, ncol), header=columnnames(cols))
        rows = Vector{Int}[Int[]]
        return cells, rows
    end

    col = cols[1]
    refs = refarray(col)
    pool = refpool(col)
    labeled = pool !== nothing && eltype(refs) <: Unsigned
    if !labeled
        refs, invpool, pool = _label(col)
    end
    mult = length(pool)
    if ncol > 1
        # Make a copy to be used as cache
        labeled && (refs = collect(refs))
        @inbounds for n in 2:ncol
            col = cols[n]
            refsn = refarray(col)
            pool = refpool(col)
            if pool === nothing || !(eltype(refsn) <: Unsigned)
                refsn, invpool, pool = _label(col)
            end
            multn = length(pool)
            _mult!(refs, mult, refsn)
            mult = mult * multn
        end
    end
    cells, rows = _cellrows(cols, _groupfind(refs))
    return cells, rows
end

findcell(names, data, esample=Colon()) =
    findcell(subcolumns(data, names, esample))
