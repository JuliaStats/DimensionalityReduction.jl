DimensionalityReduction.jl
==========================

# Algorithms

* Principal Component Analysis (PCA)

# PCA Usage

    using DimensionalityReduction

    # simulate 100 random observations
    # rotate and scale as well
    X = randn(100,2) * [0.8 0.7; 0.9 0.5]
    Xpca = pca(X)

Rows of `X` each represent a data point (i.e., a different repetition of the experiment),
and columns of `X` represent the different variables measured.

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

Centering is done by subtracting the mean, and scaling by normalizing each variable by its
standard deviation.

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

ICA has been deprecated.

# t-SNE Usage

t-SNE has been deprecated.

# NMF

NMF has been moved into a separate [package](https://github.com/JuliaStats/NMF.jl).
