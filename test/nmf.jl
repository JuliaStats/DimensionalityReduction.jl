using DimensionalityReduction

X = hcat(eye(2), eye(2))
results = nmf(X, 2)
R = results.W * results.H
@assert maximum(X - R) < 10e-4

X = hcat(eye(2), eye(2))
X = vcat(X, X, X, X)
results = nmf(X, 2)
R = results.W * results.H
@assert maximum(X - R) < 10e-4
