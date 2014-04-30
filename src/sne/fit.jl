function sne(
    X::Matrix{Float64};
    k::Integer = 2,
    iterations::Integer = 1000,
    perplexity::Real = 0.5,
    verbose::Bool = true,
    tolerance::Real = 10e-12,
)
    d, n = size(X)

    Y = 10e-4 * randn(k, n)

    D2_X = pairwise(SqEuclidean(), X)
    D2_Y = Array(Float64, n, n)

    P = Array(Float64, n, n)
    σ = Array(Float64, n)
    P_sne!(P, D2_X, σ, perplexity)

    Q = Array(Float64, n, n)

    gr = Array(Float64, k, n)
    sumsqgr = zeros(Float64, k, n)

    prevcost, cost = NaN, NaN
    for iteration in 1:iterations
        pairwise!(D2_Y, SqEuclidean(), Y)
        Q_sne!(Q, D2_Y)

        prevcost, cost = cost, f_sne(P, Q)
        if abs(cost - prevcost) < tolerance
            break
        end
        if verbose
            @printf("%d\t%.8f\t%.8f\n", iteration, cost, abs(cost - prevcost))
        end
        g_sne!(gr, P, Q, Y)

        for index in 1:(k * n)
            sumsqgr[index] += gr[index]^2
            α = 1 / sqrt(sumsqgr[index])
            Y[index] = Y[index] - α * gr[index]
        end
    end

    return Y
end

function tsne(
    X::Matrix{Float64};
    k::Integer = 2,
    iterations::Integer = 1000,
    perplexity::Real = 0.8,
    verbose::Bool = true,
    tolerance::Real = 10e-12,
)
    d, n = size(X)

    Y = 10e-4 * randn(k, n)

    D2_X = pairwise(SqEuclidean(), X)
    D2_Y = Array(Float64, n, n)

    P = Array(Float64, n, n)
    σ = Array(Float64, n)
    P_tsne!(P, D2_X, σ, perplexity)

    Q = Array(Float64, n, n)

    gr = Array(Float64, k, n)
    sumsqgr = zeros(Float64, k, n)

    prevcost, cost = NaN, NaN
    for iteration in 1:iterations
        pairwise!(D2_Y, SqEuclidean(), Y)
        Q_tsne!(Q, D2_Y)

        prevcost, cost = cost, f_sne(P, Q)
        if verbose
            @printf("%d\t%.8f\t%.8f\n", iteration, cost, abs(cost - prevcost))
        end
        if abs(cost - prevcost) < tolerance
            break
        end
        g_tsne!(gr, P, Q, Y, D2_Y)

        for index in 1:(k * n)
            sumsqgr[index] += gr[index]^2
            α = 1 / sqrt(sumsqgr[index])
            Y[index] = Y[index] - α * gr[index]
        end
    end

    return Y
end
