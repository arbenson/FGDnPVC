include("algorithms.jl")

using Base.Threads
using MAT

using ScikitLearn
@sk_import metrics: average_precision_score

"""
recovery_over_time(dataset::String, interval_in_days::Int64=10)
---------------------

Computes the precision at core size, area under the 

using the following algorithms:
- Union of minimal verex covers (UMVC)
- Degree ordering
- Betweenness centrality score ordering
- Borgatti & Everret core scores

recovery_over_time(dataset::String, interval_in_days::Int64=10)

Input parameters:
- dataset::String: dataset name
- interval_in_days::Int64: sampling interval

stores results in output/$dataset-temporal-perf-stats.mat
"""
function recovery_over_time(dataset::String, interval_in_days::Int64=10)
    srand(1234)  # for consistency
    
    TD = read_temporal_data(dataset)
    A = TemporalData2SimpleGraph(TD)
    n = size(A, 1)
    core01 = read_temporal_core(dataset, n)
    core = find(core01 .> 0)
    nc = length(core)
    
    total_days = (maximum(TD.times) - minimum(TD.times)) / 86400
    increment = interval_in_days / total_days
    
    ps = collect(increment:increment:1.0)
    np = length(ps)

    perf_deg_pacs = zeros(Float64, np)
    perf_mvc_pacs = zeros(Float64, np)
    perf_upb_pacs = zeros(Float64, np)
    perf_btw_pacs = zeros(Float64, np)
    perf_bev_pacs = zeros(Float64, np)

    shuffled_inds = shuffle(collect(1:np))
    
    all_deg_scores = zeros(Float64, length(core01), np)
    all_mvc_scores = zeros(Float64, length(core01), np)
    all_btw_scores = zeros(Float64, length(core01), np)
    all_upb_scores = zeros(Float64, length(core01), np)
    all_bev_scores = zeros(Float64, length(core01), np)

    nedges = nnz(A) / 2
    nnodes = sum(vec(sum(A, 2)) .> 0)
    frac_edges = zeros(Float64, np)
    frac_core  = zeros(Float64, np)
    frac_nodes = zeros(Float64, np)
    
    Threads.@threads for ii = 1:length(shuffled_inds)
        i = shuffled_inds[ii]
        p = ps[i]
        if Threads.threadid() == 1
            print("$(ii) of $(np)... \r")
            flush(STDOUT)
        end
        # Collect data and scores
        Ap = TemporalData2SimpleGraph(quantiled_data(TD, p), n)
        deg_order, d = degree_order(Ap)
        mvc_order = UMVC_order(Ap)
        btw_order, btw_scores = betweenness_order(Ap)
        bev_order, bev_scores = BorgattiEverett(Ap)
        non_iso_core = core[find(d[core] .> 0)]

        # precision @ core size
        end_ind = min(nc, size(Ap, 1))
        pacs(ord::Vector{Int64}) = length(intersect(ord[1:end_ind], core)) / nc
        perf_deg_pacs[i] = pacs(deg_order)
        perf_mvc_pacs[i] = pacs(mvc_order)
        perf_btw_pacs[i] = pacs(btw_order)
        perf_bev_pacs[i] = pacs(bev_order)
        perf_upb_pacs[i] = length(ni_core) / nc

        # upper bound
        upb_scores = rand(Float64, n)
        mval = maximum(upb_scores)
        upb_scores[ni_core] = mval * 10.0
        deg_scores = convert(Vector{Float64}, d)

        # Store data to avoid threading issues
        all_deg_scores[:, i] = deg_scores
        all_mvc_scores[:, i] = mvc_scores
        all_btw_scores[:, i] = btw_scores
        all_bev_scores[:, i] = bev_scores
        all_upb_scores[:, i] = upb_scores

        frac_edges[i] = (nnz(Ap) / 2) / nedges
        frac_core[i]  = length(ni_core) / nc
        frac_nodes[i] = sum(vec(sum(Ap, 2)) .> 0) / nnodes
    end

    perf_deg_auprc = zeros(Float64, np)
    perf_mvc_auprc = zeros(Float64, np)
    perf_bev_auprc = zeros(Float64, np)    
    perf_upb_auprc = zeros(Float64, np)
    perf_btw_auprc = zeros(Float64, np)
    for i = 1:np
        perf_deg_auprc[i] =
            average_precision_score(core01, all_deg_scores[:, i])
        perf_mvc_auprc[i] =
            average_precision_score(core01, all_mvc_scores[:, i])
        perf_btw_auprc[i] =
            average_precision_score(core01, all_btw_scores[:, i])
        perf_bev_auprc[i] =
            average_precision_score(core01, all_bev_scores[:, i])
        perf_upb_auprc[i] =
            average_precision_score(core01, all_upb_scores[:, i])
    end

    matwrite("output/$dataset-temporal-perf-stats.mat",
             Dict("ps"         => ps,
                  "interval"   => interval_in_days,                  
                  "frac_core"  => frac_core,
                  "frac_nodes" => frac_nodes,
                  "frac_edges" => frac_edges,
                  "ncore"      => nc,
                  "nedges"     => nedges,
                  "nnodes"     => nnodes,
                  "deg_auprc"  => perf_deg_auprc,
                  "mvc_auprc"  => perf_mvc_auprc,
                  "btw_auprc"  => perf_btw_auprc,                  
                  "upb_auprc"  => perf_upb_auprc,
                  "bev_auprc"  => perf_bev_auprc,
                  "deg_pacs"   => perf_deg_pacs,
                  "mvc_pacs"   => perf_mvc_pacs,
                  "btw_pacs"   => perf_btw_pacs,
                  "upb_pacs"   => perf_upb_pacs,
                  "bev_pacs"   => perf_bev_pacs,))
end
