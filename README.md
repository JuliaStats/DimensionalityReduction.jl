DimensionalityReduction.jl
==========================

# Algorithms

* Principal Component Analysis (PCA)
* Independant Component Analysis (ICA)
* Non-negative Matrix Factorization (NMF)
* t-Distributed Stochastic Neighbor Embedding (t-SNE)

# PCA Usage

    using DimensionalityReduction

    # simulate 100 random observations
    # rotate and scale as well
    X = randn(100,2) * [0.8 0.7; 0.9 0.5]
    Xpca = pca(X)

Attributes:

    Xpca.rotation                # principal components
    Xpca.scores                  # rotated X
    Xpca.standard_deviations     # square roots of the eigenvalues
    Xpca.proportion_of_variance  # fraction of variance brought by each principal component
    Xpca.cumulative_variance     # cumulative proportion of variance

By default, `pca()` uses SVD decomposition. Alternatively, `pcaeig(X)` will calculate
directly the eigenvectors of the covariance matrix.

`pca()` centers and re-scales input data by default.
This is controlled by the `center` and `scale` keyword arguments:

	pca(X::Matrix ; center::Bool, scale::Bool)

If `scale` is true (default), then the principal components of the data are also
scaled back to the original space and saved to `Xpca.rotation`

To overlay the principal components on top of the data with [PyPlot](https://github.com/stevengj/PyPlot.jl)

	using PyPlot
	plot( X[:,1], X[:,2], "r." )  # point cloud

	# get data center
	ctr = mean( X, 1 )

	# plot principal components as lines
	#  weight by their standard deviation
	PCs = Xpca.rotation
	for v=1:2
		weight = Xpca.standard_deviations[v]
		plot( ctr[1] + weight * [0, PCs[1,v]], 
			  ctr[2] + weight * [0, PCs[2,v]],
			  linewidth = 2)
	end



To make a biplot with [PyPlot](https://github.com/stevengj/PyPlot.jl)

	using PyPlot
	scores = Xpca.scores[:,1:2]
	plot( scores[:,1], scores[:,2], "r." )


To make a biplot with [Gadfly](http://dcjones.github.io/Gadfly.jl/):

    using Gadfly
    scores = Xpca.scores[:,1:2]
    pl = plot(x=scores[:,1],y=scores[:,2], Geom.point)
    draw(PNG("pca.png", 6inch, 6inch), pl)

Starting from a DataFrame:

    using RDatasets
    iris = data("datasets", "iris")
    iris = convert(Array,DataArray(iris[:,1:4]))
    Xpca = pca(iris)

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
