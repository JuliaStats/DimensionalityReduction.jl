DimensionalityReduction.jl
==========================

# Algorithms

* PCA
* ICA
* NMF

# PCA Usage

	using DimensionalityReduction

	srand(1)

	X = hcat(10 * randn(10), randn(10))
	cov(X)

	theta = pi / 4.0
	R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
	X = X * R
	cov(X)

	results = pca(X)

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
