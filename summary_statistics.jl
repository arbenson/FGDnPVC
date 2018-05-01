include("common.jl")

using MathProgBase
using Gurobi

"""
minimum_VC_size
---------------

Use integer programming formulation to find the minimum vertex cover size.

size = minimum_VC_size(A::SpIntMat, upper_bound::Int64)

Input parameters:
- A::SpIntMat: adjacency matrix of simple undirected graph
- upper_bound::Int64: upper bound on the minimum vertex cover size
"""
function minimum_VC_size(A::SpIntMat, upper_bound::Int64)
    n = size(A, 1)
    vartypes = Symbol[]
    for _ in 1:n; push!(vartypes, :Bin); end
    T = triu(A, 1)
    I = Int64[]
    J = Int64[]
    edge_ind = 1
    row_lb = Int64[]
    row_ub = Int64[]
    for i in 1:size(T, 2), j in nz_row_inds(T, i)
        if T[j, i] > 0
            push!(I, edge_ind, edge_ind)
            push!(J, i, j)
            edge_ind += 1
            push!(row_ub, 2)
            push!(row_lb, 1)
        end
    end
    c = ones(Int64, n)
    # append one for the upper bound
    append!(J, collect(1:n))
    append!(I, ones(Int64, n) * edge_ind)
    push!(row_lb, 1)
    push!(row_ub, upper_bound)
    B = convert(SpIntMat, sparse(I, J, ones(length(I)), edge_ind, n))
    lb = 0
    ub = 1
    sol = mixintprog(c, B, row_lb, row_ub, vartypes, lb, ub, GurobiSolver())
    return convert(Int64, sol.objval)
end

"""
summary_statistics
------------------

Compute summary statistics of a dataset.

summary_statistics(dataset::String)

Input parameters:
- dataset::String: name of dataset
"""
function summary_statistics(dataset::String)
    TD =  read_temporal_data(dataset)
    total_days = convert(Int64, round((maximum(TD.times) - minimum(TD.times)) / 86400))
    A = TemporalData2SimpleGraph(TD)
    degs = vec(sum(A, 2))
    n = size(A, 1)
    nnodes = sum(degs .> 0)
    core01 = read_core(dataset, n)
    core = find(core01 .== 1)
    k = sum(degs[core] .> 0)
    kstar = minimum_VC_size(A, k)
    kstar = minimum_VC_size(A::SpIntMat, k)

    # Find all nodes that connect outside the core
    core_connecting = zeros(Int64, n)
    for c in core
        for j in nz_row_inds(A, c)
            if core01[j] == 0
                core_connecting[c] = 1
                break
            end
        end
    end

    in_interior = zeros(Int64, n)    
    for c in core
        if core_connecting[c] == 0
            # Check if node is in the interior
            s = 0
            for j in nz_row_inds(A, c)
                assert(core01[j] == 1)
                if core_connecting[j] == 1
                    s += 1
                    break
                end
            end
            if s == 0
                in_interior[c] = 1
            end
        end
    end

    interior_connecting = zeros(Int64, n)
    for c in core
        if in_interior[c] == 1
            for j in nz_row_inds(A, c)
                interior_connecting[j] = 1
            end
            #interior_connecting[c] = 1
        end
    end
    nedges = nnz(A) / 2

    println("dataset & n & m & time span (days) & |C| & k* & Bound 1 & Bound 2 & fraction C connecting to periphery & fraction of C with 2-hop nbrhood in C")
    println(@sprintf("%s & %d & %d & %d & %d & %d & %0.3f & %0.3f & %0.2f & %0.2f",
                     dataset,
                     nnodes,
                     nedges,
                     total_days,
                     k,
                     kstar,
                     ((k + 1)^2 / 4 + k) / nnodes,
                     (k - kstar + 2) * kstar / nnodes,
                     sum(core_connecting) / length(core),
                     sum(interior_connecting) / length(core) ))
end
