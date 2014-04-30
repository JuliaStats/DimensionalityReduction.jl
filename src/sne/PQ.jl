function P_sne!(
    P::Matrix{Float64},
    D2_X::Matrix{Float64},
    σ::Vector{Float64},
    perplexity::Real,
)
    n, n = size(D2_X)
    target = perplexity * 1.0 + (1.0 - perplexity) * (n - 1.0)
    for j in 1:n
        lower, upper = 10e-16, 10e16
        σ_j = (lower + upper) / 2
        for pass in 1:100
            normalizer = 0.0
            for i in 1:n
                if i != j
                    P_ij = exp(-D2_X[i, j] / σ_j)
                    P[i, j] = P_ij
                    normalizer += P_ij
                else
                    P[i, j] = 0.0
                end
            end
            entropy = 0.0
            for i in 1:n
                if i != j
                    P_ij = P[i, j] / normalizer
                    P[i, j] = P_ij
                    entropy -= P_ij * log2(P_ij)
                end
            end
            perplexity = 2^entropy
            if isapprox(perplexity, target)
                break
            else
                if target < perplexity
                    upper = σ_j
                    σ_j = (lower + upper) / 2
                else
                    lower = σ_j
                    σ_j = (lower + upper) / 2
                end
            end
        end
        σ[j] = σ_j
    end
    return
end

function P_tsne!(
    P::Matrix{Float64},
    D2_X::Matrix{Float64},
    σ::Vector{Float64},
    perplexity::Real,
)
    P_sne!(P, D2_X, σ, perplexity)
    n, n = size(P)
    for j in 1:n
        for i in (j + 1):n
            P_ij = P[i, j]
            P_ji = P[j, i]
            tmp = (P_ij + P_ji) / (2 * n)
            P[i, j] = tmp
            P[j, i] = tmp
        end
    end
    return
end

function Q_sne!(
    Q::Matrix{Float64},
    D2_Y::Matrix{Float64},
)
    n, n = size(D2_Y)
    for j in 1:n
        normalizer = 0.0
        for i in 1:n
            if i != j
                Q_ij = exp(-D2_Y[i, j])
                Q[i, j] = Q_ij
                normalizer += Q_ij
            else
                Q[i, j] = 0.0
            end
        end
        for i in 1:n
            if i != j
                Q[i, j] = Q[i, j] / normalizer
            end
        end
    end
    return
end

function Q_tsne!(
    Q::Matrix{Float64},
    D2_Y::Matrix{Float64},
)
    n, n = size(D2_Y)
    normalizer = 0.0
    for j in 1:n
        for i in 1:n
            if i != j
                Q_ij = 1 / (1 + D2_Y[i, j])
                Q[i, j] = Q_ij
                normalizer += Q_ij
            else
                Q[i, j] = 0.0
            end
        end
    end
    for j in 1:n
        for i in 1:n
            if i != j
                Q[i, j] = Q[i, j] / normalizer
            end
        end
    end
    return
end
