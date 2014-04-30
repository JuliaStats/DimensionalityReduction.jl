module TestSNE
    using DimensionalityReduction

    X = hcat(randn(10, 100), randn(10, 100) .+ 10)

    Y = sne(X, verbose = false)
    Y = tsne(X, verbose = false)
end
