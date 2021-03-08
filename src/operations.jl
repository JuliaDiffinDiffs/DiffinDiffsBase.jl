# Obtain unique labels for all row-wise pairs of values from a1 and a2
function _relabel!(a1::Array, mult1::Integer, a2::AbstractArray)
    size(a1) == size(a2) || throw(DimensionMismatch(
        "cannot match array of size $(size(a1)) with array of size $(size(a2))"))
    refs, invpool, pool = _label(a2)
    _mult!(a1, mult1, refs)
    return mult1 * length(pool)
end

function _relabel!(a1::Array, mult1::Integer, a2::PooledArray)
    size(a1) == size(a2) || throw(DimensionMismatch(
        "cannot match array of size $(size(a1)) with array of size $(size(a2))"))
    _mult!(a1, mult1, a2.refs)
    return mult1 * length(a2.pool)
end

function _mult!(a1::Array, mult1::Integer, a2::Array)
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

# Convert the keys from refs to cells
function _cellrows(cols::SubColumns, refrows::IdDict)
    ncol = size(cols, 2)
    cellrows = IdDict{Tuple, Vector{Int}}()
    for rows in values(refrows)
        cell = ntuple(n->cols[n][rows[1]], ncol)
        cellrows[cell] = rows
    end
    return cellrows
end

"""
    findcell(cellnames, data, esample=Colon())

Group the row indices of a subsample of `data` over `esample`
so that the row-wise combinations of values from columns indexed by `cellnames`
are the same within each group.

Note that unless `esample` covers all rows of `data`,
the row indices are those for the subsample selected based on `esample`
rather than those for the full `data`.

# Returns
- `IdDict{Tuple, Vector{Int}}`: a map from row-wise combinations of values to row indices of these combinations.
"""
function findcell(cellnames, data, esample=Colon())
    cols = SubColumns(data, cellnames, esample)
    isempty(cols) && throw(ArgumentError("empty data columns"))
    ncol = size(cols, 2)
    pooled = cols[1] isa PooledArray
    if pooled
        refs = cols[1].refs
        mult = length(cols[1].pool)
    else
        refs, invpool, pool = _label(cols[1])
        mult = length(pool)
    end
    if ncol > 1
        pooled && (refs = copy(refs))
        for n in 2:ncol
            mult = _relabel!(refs, mult, cols[n])
        end
    end
    cellrows = _cellrows(cols, _groupfind(refs))
    return cellrows
end
