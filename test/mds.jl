using DimensionalityReduction

X = eye(4)
X[1, 1] = 2.0

results = mds(dist(X), 2)
@assert size(results.X) == (4, 2)
