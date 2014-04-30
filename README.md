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

# SNE/t-SNE Usage

SNE and t-SNE are two closely related algorithms that can be used to produce
two-dimensional visualizations of high-dimensional datasets. They are available
through the `sne` and `tsne` functions:

```
X = hcat(randn(10, 100), randn(10, 100) .+ 10)

Y = sne(X)
Y = tsne(X)
```

Both algorithms expose a variety of options. These are made available using
keyword arguments:

* `k` is the number of dimensions in which the output should lie. It defaults
   to `2`
* `iterations` determines the number of iterations of gradient descent used
   when fitting the model. It defaults to `1000`, which can be prohibitive
   for large data sets.
* `perplexity` determines whether the model focused on local or global
   structure when embedding. It can take on any value strictly between `0` and
   `1` and defaults to `0.9`.
* `verbose` determines whether the gradient descent process will print
   out verbose information about the cost function value's at each iteration.
   It defaults to `true` because the algorithm is slow and may seem to hang.
* `tolerance` determines the minimum change in the cost function before
   the gradient descent algorithm will be terminated. It defaults to `10e-12`.

# NMF

NMF has been moved into a separate [package](https://github.com/JuliaStats/NMF.jl).
