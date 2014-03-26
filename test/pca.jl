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


###  Tests related to issue #13  ###

# to compare columns of A and B, regardless of sign swapping
function mycomp( A, B ) 
  err = 0.0
  for i=1:size(A,1)
    err += min( sum(abs( A[:,i] - B[:,i] )), sum(abs( A[:,i] + B[:,i] )) )
  end

  return err
end

srand(1238)
X = randn(200,4) * rand(4,4)

testMatrix = { mycomp( pcasvd(X, scale=sc, center=ctr).rotation, 
				pcaeig(X, scale=sc, center=ctr).rotation )
  for sc=[true,false], ctr=[true,false] }

@assert all( testMatrix .< 1e-12 )
