DimensionalityReduction.jl
==========================

# Algorithms

* PCA
* ICA
* NMF

# Usage

	require("DimensionalityReduction")
	using DimensionalityReduction

	srand(1)

	x = hcat(10 * randn(10), randn(10))
	cov(x)

	# Rotate through an angle to increase covariance
	theta = pi / 4.0
	r = [cos(theta) -sin(theta); sin(theta) cos(theta)]
	x = (r * x')'
	cov(x)

	results = pca(x)
