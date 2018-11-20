using LinearAlgebra
using SparseArrays

const SpIntMat = SparseMatrixCSC{Int64,Int64}
nz_row_inds(A::SpIntMat, ind::Int64) = A.rowval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]
nz_row_vals(A::SpIntMat, ind::Int64) = A.nzval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]

struct TemporalData
    sources::Vector{Int64}
    dests::Vector{Int64}
    times::Vector{Float64}
end

"""
TemporalData2SimpleGraph
------------------------

Convert a temporal graph dataset to a simple, undirected, unweighted graph.

A = TemporalData2SimpleGraph(data::TemporalData,
                             n::Int64=max(maximum(data.sources),
                                          maximum(data.dests)))

Input parameters:
- data::TemporalData: temporal graph dataset
- n::Int64: number of nodes

returns a sparse adjacency matrix
"""
function TemporalData2SimpleGraph(data::TemporalData,
                                  n::Int64=max(maximum(data.sources),
                                               maximum(data.dests)))
    I = data.sources
    J = data.dests
    A = convert(SpIntMat, sparse(I, J, ones(length(I)), n, n))
    A = max.(A, A')
    A = min.(A, 1)
    A -= Diagonal(A)
    return A
end

"""
quantiled_data
--------------

Get first p fraction of temporal data.

TD = quantiled_data(data::TemporalData, p::Float64)

Input parameters:
- data::TemporalData: temporal graph dataset
- p::Float64: fraction of data to keepxsxs

returns a new TemporalData object
"""
function quantiled_data(data::TemporalData, p::Float64)
    new_sources = Int64[]
    new_dests   = Int64[]
    new_times   = Float64[]
    q = quantile(data.times, p)
    for (s, d, t) in zip(data.sources, data.dests, data.times)
        if t <= q
            push!(new_sources, s)
            push!(new_dests,   d)
            push!(new_times,   t)
        end
    end
    return TemporalData(new_sources, new_dests, new_times)
end

"""
read_core
---------

Read the core of a dataset.

core01 = read_core(dataset::String, num_nodes::Int64)

Input parameters:
- dataset::String: dataset name
- num_nodes::Int64: number of nodes in the graph

returns a length-num_nodes vector which has value 1 in index i if node i is in
the core and value 0 if node i is not in the core.
"""
function read_core(dataset::String, num_nodes::Int64)
    core = zeros(Int64, num_nodes)
    core_ids = convert(Vector{Int64},
                       readdlm("data/$dataset/core-$dataset.txt")[:, 1])
    core[core_ids] .= 1
    return core
end

"""
read_temporal_data
------------------

Read a temporal graph from a text file, where each row is

node node timestamp


TD = read_temporal_data(dataset::String)

Input parameters:
- dataset::String: dataset name

returns a TemporalData object
"""
function read_temporal_data(dataset::String)
    sources = Int64[]
    dests   = Int64[]
    times   = Float64[]
    open("data/$dataset/$dataset.txt") do f
        for line in eachline(f)
            data = split(line)
            push!(sources, parse(Int64,   data[1]))
            push!(dests,   parse(Int64,   data[2]))
            push!(times,   parse(Float64, data[3]))
        end
    end
    return TemporalData(sources, dests, times)
end
