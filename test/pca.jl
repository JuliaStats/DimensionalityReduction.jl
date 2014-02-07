using DimensionalityReduction

X = hcat(eye(2), eye(2))
results = pcaeig(X ; center=false, scale=false )
@assert maximum(X - results.scores * results.rotation') < 10e-4
results = pcasvd(X ; center=false, scale=false)
@assert maximum(X - results.scores * results.rotation') < 10e-4
results = pca(X ; center=false, scale=false)
@assert maximum(X - results.scores * results.rotation') < 10e-4

X = hcat(eye(2), eye(2))
X = vcat(X, X, X, X)
results = pcaeig(X ; center=false, scale=false)
@assert maximum(X - results.scores * results.rotation') < 10e-4
results = pcasvd(X ; center=false, scale=false)
@assert maximum(X - results.scores * results.rotation') < 10e-4
results = pca(X ; center=false, scale=false)
@assert maximum(X - results.scores * results.rotation') < 10e-4

srand(1)

X = hcat(10 * randn(10), randn(10))
cov(X)

theta = pi / 4.0
R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
X = X * R
cov(X)

results = pca(X)
