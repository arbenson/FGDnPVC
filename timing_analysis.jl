include("algorithms.jl")

"""
timing_experiment
-----------------

Script to run timing experiments.

timing_experiment(dataset::String, alg::String)

Input parameters:
- dataset::String: name of dataset
- alg::String: which algorithm to time 
               ("degree", "UMVC", "betweenness", or "BorgattiEverett")
"""
function timing_experiment(dataset::String, alg::String)
    A = TemporalData2SimpleGraph(read_temporal_data(dataset))
    if     alg == "degree";          @time degree_order(A)
    elseif alg == "UMVC";            @time UMVC_order(A)
    elseif alg == "betweenness";     @time betweenness_order(A)
    elseif alg == "BorgattiEverett"; @time BorgattiEverett(A)
    else throw("Unknown algorithm $alg"); end
end
