include("common.jl")
using LightGraphs

function UMVC(A::SpIntMat, num_iters::Int64=300)
    edge_vec = [(i, j) for (i, j) in zip(findnz(A)[1:2]...) if i < j]
    n = size(A, 1)
    umvc_vec = zeros(Int64, n)
    for _ in 1:num_iters
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
            for c in shuffle(find(cover .== 1))
                nbrs = nz_row_inds(A, c)
                if sum(cover[nbrs]) == length(nbrs)
                    cover[c] = 0
                    reduced = true
                end
            end
            if !reduced; break; end
        end
        umvc_vec[cover .== 1] = 1
    end

    # Order and score the nodes
    d = vec(sum(A, 2))  # degrees of nodes
    return [sort(find(umvc_vec .== 1), by=v->d[v], rev=true);
            sort(find(umvc_vec .== 0), by=v->d[v], rev=true)]
end

function degree_order(A::SpIntMat)
    d = vec(sum(A, 2))
    n = size(A, 1)
    return (sort(collect(1:n), by=v->d[v], rev=true), d)
end

function BorgattiEverett(A::SpIntMat, max_iter::Int64=10000, tol::Float64=1e-10)
    n = size(A, 1)
    d = vec(sum(A, 2))
    c = rand(n)
    c[d .== 0] = 0
    c /= norm(c, 2)
    for iter = 1:max_iter
        num = A * c
        denom = sum(c .^ 2) - c .^ 2
        next_c = num ./ denom
        next_c /= norm(next_c, 2)
        diff = norm(next_c - c, 2)
        c = next_c
        if diff < tol; break; end
    end
    return (sort(collect(1:n), by= v -> c[v], rev=true), c)
end

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
