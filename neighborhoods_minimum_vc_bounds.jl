include("common.jl")

using Base.Threads
using MAT

"""
read_simple_graph_txt
---------------------

Read a simple, undirected, unweighted graph from a text file. 
Example datasets are in the directory data/static/

(A, index_map) = read_simple_graph_txt(dataset::String)

Input parameters:
- dataset::String: name of dataset
"""
function read_simple_graph_txt(dataset::String)
    filename = "data/static/$dataset.txt"
    index_map = Dict{Int64,Int64}()
    index_map_vec = Int64[]
    function get_mapped_index(x::Int64)
        if !haskey(index_map, x)
            next_index = length(index_map) + 1
            index_map[x] = next_index
            push!(index_map_vec, x)
            return next_index
        end
        return index_map[x]
    end

    # Read data
    I = Int64[]
    J = Int64[]
    open(filename) do f
        for line in eachline(f)
            # Skip lines starting with '#' or '%'
            if line[1] == '#' || line[1] == '%'; continue; end
            edge = split(line)
            u = parse(Int64, edge[1])
            v = parse(Int64, edge[2])
            push!(I, get_mapped_index(u))
            push!(J, get_mapped_index(v))
        end
    end

    # Form adjacency matrix
    n = max(maximum(I), maximum(J))
    A = convert(SpIntMat, sparse(I, J, ones(length(I)), n, n))
    A = max.(A, A')
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    return (A, index_map_vec)
end

"""
matching_approx_bounds
----------------------

Get approximation bounds from num_iters random initializations of the maximal
matching 2-approximation for the minimum vertex cover problem on the undirected
graph with adjacency matrix A.

(u, l) = matching_approx_bounds(SpIntMat, num_iters::Int64=20)

Input parameters:
- SpIntMat: the adjacency matrix of an undirected graph
- num_iters::Int64=20: number of random iterations of the approximation

returns tuple of integers (u, l), which are the worst upper bound and best lower
bound on the minimum vertex cover size.
"""
function matching_approx_bounds(A::SpIntMat, num_iters::Int64=20)
    edge_vec = [(i, j) for (i, j) in zip(findnz(A)[1:2]...) if i < j]
    n = size(A, 1)
    upper_bounds = Int64[]
    lower_bounds = Int64[]
    for _ in 1:num_iters
        # Run 2-approximation
        cover = zeros(Int64, n)
        edge_queue = copy(shuffle(edge_vec))
        while length(edge_queue) > 0
            i, j = pop!(edge_queue)
            if cover[[i, j]] == [0, 0]; cover[[i, j]] = [1, 1]; end
        end
        push!(upper_bounds, sum(cover) / 2)  # bad bound
        # Reduce to a minimal cover
        while true
            reduced = false
            for c in shuffle(find(cover .== 1))
                nbrs = nz_row_inds(A, c)
                if sum(cover[nbrs]) == length(nbrs)
                    cover[c] = 0
                    reduced = true
                end
            end
            if !reduced; break; end
        end
        push!(lower_bounds, sum(cover))  # good bound
    end
    return (maximum(upper_bounds), minimum(lower_bounds))
end

"""
neighborhoods_minimum_VC_bounds
-------------------------------

Collect bounds on the minimum vertex cover size of 1-hop neighborhood covers of
2-hop neighborhoods (ego excluded).

neighborhoods_minimum_VC_bounds(dataset::String, niter::Int64=20)

Input parameters:
- dataset::String: name of dataset
- niter::Int64=20: number of interations of the maximal matching algorithm to use

Writes output file output/$dataset-neighborhood-stats.mat.
"""
function neighborhoods_minimum_VC_bounds(dataset::String, niter::Int64=20)
    A = read_simple_graph_txt(dataset)[1]

    n = size(A, 2)
    degs = vec(sum(A, 2))
    N1_sizes = copy(degs)
    N2_sizes = zeros(Int64, n)
    counted  = zeros(Int64, n)
    uvals   = zeros(Int64, n)
    lvals   = zeros(Int64, n)    

    Threads.@threads for i = 1:n
        if Threads.threadid() == 1
            print("$(i) of $n \r")
            flush(STDOUT)
        end
        nbrs = nz_row_inds(A, i)
        k = length(nbrs)
        if (k + 1)^2 / 4 + k <= n  # only include if bound is nontrivial
            counted[i] = 1
            all = copy(nbrs)
            for nbr in nbrs
                append!(all, [u for u in nz_row_inds(A, nbr) if u != i])
            end
            N2_sizes[i] = length(unique(all))

            # Form adjacency matrix for 2-hop subgraph
            I = Int64[]
            J = Int64[]
            for nbr in nbrs, j in nz_row_inds(A, nbr)
                if j != i
                    push!(I, min(nbr, j))
                    push!(J, max(nbr, j))
                end
            end
            if length(I) == 0
                # If there are no edges, then nothing is done.nn
                uvals[i] = length(nbrs)
                lvals[i] = length(nbrs)                
            else
                # Run the 
                bn = max(maximum(I), maximum(J))
                B = convert(SpIntMat, sparse(I, J, ones(size(I)), bn, bn))
                B = max.(B, B')
                approx = matching_approx_bounds(B, niter)
                uvals[i] = approx[1]
                lvals[i] = approx[2]
            end
        end
    end

    matwrite("output/$dataset-neighborhood-stats.mat",
             Dict("counted" => counted,
                  "N1"      => N1_sizes,
                  "N2"      => N2_sizes,
                  "u"       => uvals,
                  "l"       => lvals))
end
