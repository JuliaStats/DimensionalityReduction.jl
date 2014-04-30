function f_sne(
    P::Matrix{Float64},
    Q::Matrix{Float64},
)
    n, n = size(P)
    kl = 0.0
    for j in 1:n
        for i in 1:n
            if i != j
                kl += P[i, j] * log(P[i, j] / Q[i, j])
            end
        end
    end
    return kl
end

function f_tsne(
    P::Matrix{Float64},
    Q::Matrix{Float64},
)
    n, n = size(P)
    kl = 0.0
    for j in 1:n
        for i in 1:n
            if i != j
                kl += P[i, j] * log(P[i, j] / Q[i, j])
            end
        end
    end
    return kl
end

function g_sne!(
    gr::Matrix{Float64},
    P::Matrix{Float64},
    Q::Matrix{Float64},
    Y::Matrix{Float64},
)
    k, n = size(Y)
    fill!(gr, 0.0)
    for dim in 1:k
        for i in 1:n
            for j in 1:n
                if i != j
                    prdiff = (P[j, i] - Q[j, i] + P[i, j] - Q[i, j])
                    gr[dim, i] += 2 * prdiff * (Y[dim, i] - Y[dim, j])
                end
            end
        end
    end
    return
end

function g_tsne!(
    gr::Matrix{Float64},
    P::Matrix{Float64},
    Q::Matrix{Float64},
    Y::Matrix{Float64},
    D2_Y::Matrix{Float64},
)
    k, n = size(Y)
    fill!(gr, 0.0)
    for dim in 1:k
        for i in 1:n
            for j in 1:n
                if i != j
                    prdiff = (P[i, j] - Q[i, j])
                    tkern = 1 / (1 + D2_Y[i, j])
                    gr[dim, i] += 4 * prdiff * (Y[dim, i] - Y[dim, j]) * tkern
                end
            end
        end
    end
    return
end
