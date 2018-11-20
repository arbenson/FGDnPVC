include("common.jl")

using LightGraphs  # for betweenness centrality

"""
UMVC_order
----------

Uses the union of minimal vertex cover algorithm to compute an ordering
of nodes in terms of favorability of being in the core.

order = UMVC_order(A::SpIntMat, ncovers::Int64=300)

Input parameters:
- A::SpIntMat: adjacency matrix
- ncovers::Int64=300: number of minimal vertex covers to use

returns an ordering of the nodes in decreasing favorability
of including in the core
"""
function UMVC_order(A::SpIntMat, ncovers::Int64=300)
    edge_vec = [(i, j) for (i, j) in zip(findnz(A)[1:2]...) if i < j]
    n = size(A, 1)
    umvc_vec = zeros(Int64, n)
    for _ in 1:ncovers
        # Run 2-approximation
        cover = zeros(Int64, n)
        edge_queue = copy(shuffle(edge_vec))
        while length(edge_queue) > 0
            i, j = pop!(edge_queue)
            if cover[[i, j]] == [0, 0]; cover[[i, j]] = [1, 1]; end
        end
        # Reduce to a minimal cover
        while true
            reduced = false
            for c in shuffle(findall(cover .== 1))
                nbrs = nz_row_inds(A, c)
                if sum(cover[nbrs]) == length(nbrs)
                    cover[c] = 0
                    reduced = true
                end
            end
            if !reduced; break; end
        end
        umvc_vec[cover .== 1] .= 1
    end

    # Order and score the nodes
    d = vec(sum(A, dims=2))  # degrees of nodes
    return [sort(findall(umvc_vec .== 1), by=v->d[v], rev=true);
            sort(findall(umvc_vec .== 0), by=v->d[v], rev=true)]
end

"""
degree_order
------------

Computes the ordering of nodes in terms of their degree.

(order, d) = degree_order(A::SpIntMat)

Input parameters:
- A::SpIntMat: adjacency matrix

returns an (order, degree) tuple of the order of the nodes by decreasing
degree along with the actual degrees.
"""
function degree_order(A::SpIntMat)
    d = vec(sum(A, dims=2))
    n = size(A, 1)
    return (sort(collect(1:n), by=v->d[v], rev=true), d)
end

"""
betweenness_order(A::SpIntMat)
------------

Computes the ordering of nodes in terms of their betweenness centrality.

(order, btw) = betweenness_order(A::SpIntMat)

Input parameters:
- A::SpIntMat: adjacency matrix

returns an (orderg, scores) tuple of the order of the nodes by decreasing
betweenness centrality score along with the actual centrality scores.
"""
function betweenness_order(A::SpIntMat)
    G = SimpleGraph()
    for i = 1:size(A, 1); add_vertex!(G); end
    T = triu(A)
    for i in 1:size(T, 2), j in nz_row_inds(T, i)
        add_edge!(G, (i, j))
    end
    bw = betweenness_centrality(G)
    n = size(A, 1)
    return (sort(collect(1:n), by= v->bw[v], rev=true), bw)
end

"""
BorgattiEverett_order
---------------------

Computes the ordering of nodes in terms of the score proposed in
"Models of core/periphery structures", Borgatti & Everett, Social Networks, 2000.
Uses the iterative algorithm from
"The minimum residual method of factor analysis", Comrey, Psychological Reports, 1962.

(order, core_scores) = BorgattiEverett_order(A::SpIntMat, max_iter::Int64=10000, tol::Float64=1e-10)

Input parameters:
- A::SpIntMat: adjacency matrix
- max_iter::Int64: maximum number of iterations
- tol::Float64: stopping tolerance of algorithm

returns an (order, scores) tuple of the order of the nodes by decreasing
betweenness centrality score along with the actual core scores.
"""
function BorgattiEverett_order(A::SpIntMat, max_iter::Int64=10000,
                               tol::Float64=1e-10)
    n = size(A, 1)
    d = vec(sum(A, dims=2))
    c = rand(n)
    c[d .== 0] .= 0
    c /= norm(c, 2)
    for iter = 1:max_iter
        num = A * c
        denom = sum(c .^ 2) .- c .^ 2
        next_c = num ./ denom
        next_c /= norm(next_c, 2)
        diff = norm(next_c - c, 2)
        c = next_c
        if diff < tol; break; end
    end
    return (sort(collect(1:n), by= v -> c[v], rev=true), c)
end
