DimensionalityReduction.jl
==========================

# Algorithms

* Principal Component Analysis (PCA)
* Independant Component Analysis (ICA)
* Non-negative Matrix Factorization (NMF)
* t-Distributed Stochastic Neighbor Embedding (t-SNE)

# PCA Usage

    using DimensionalityReduction

    X = reshape(rand(40),10,4)
    Xpca = pca(X)

Attributes:

    Xpca.rotation                # principal components
    Xpca.scores                  # rotated X
    Xpca.standard_deviations     # square roots of the eigenvalues
    Xpca.proportion_of_variance  # fraction of variance brought by each principal component
    Xpca.cumulative_variance     # cumulative proportion of variance

By default, pca() uses SVD decomposition. Alternatively, `pcaeig(X)` will calculate
directly the eigenvectors of the covariance matrix.

# ICA Usage

    using DimensionalityReduction

    # Generate true sources
    S_true = rand(5,1000)

    # Mixing matrix
    H_true = randn(5, 5)

    # generate observed signals
    X = H_true*S_true

    results = ica(X)

# NMF Usage

    using DimensionalityReduction

    X = hcat(eye(2), eye(2))
    X = vcat(X, X, X, X)
    results = nmf(X, 2)

# t-SNE Usage

    using DimensionalityReduction

    X = hcat(eye(2), eye(2))
    X = vcat(X, X, X, X)
    results = tsne(X, 2)
