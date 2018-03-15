const SpIntMat = SparseMatrixCSC{Int64,Int64}
nz_row_inds(A::SpIntMat, ind::Int64) = A.rowval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]
nz_row_vals(A::SpIntMat, ind::Int64) = A.nzval[A.colptr[ind]:(A.colptr[ind + 1] - 1)]

immutable TemporalData
    sources::Vector{Int64}
    dests::Vector{Int64}
    times::Vector{Float64}
end

function TemporalData2Graph(data::TemporalData,
                            n::Int64=max(maximum(data.sources),
                                         maximum(data.dests)))
    I = data.sources
    J = data.dests
    return convert(SpIntMat, sparse(I, J, ones(length(I)), n, n))
end

function TemporalData2SimpleGraph(data::TemporalData,
                                  n::Int64=max(maximum(data.sources),
                                               maximum(data.dests)))
    A = TemporalData2Graph(data, n)
    A = max.(A, A')
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    return A
end

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

function read_core(dataset::String, num_nodes::Int64)
    core = zeros(Int64, num_nodes)
    core_ids = convert(Vector{Int64},
                       readdlm("data/$dataset/core-$dataset.txt")[:, 1])
    core[core_ids] = 1
    return core
end

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
