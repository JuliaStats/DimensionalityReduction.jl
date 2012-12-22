function eigenvector_pca{T <: Real}(x::Matrix{T})
  c_x = cov(x)
  (l, z) = eig(c_x)
  p = z'
  # eigenvectors are in reverse order
  l = reverse(l)
  p = p[:, reverse(1:size(p, 2))]
  y = x * p
  return PrincipalComponentsDecomposition(p, y, sqrt(l), l / sum(l), cumsum(l) / sum(l))
end

function svd_pca{T <: Real}(x::Matrix{T})
  (a, b, c) = svd(cov(x))
  return PrincipalComponentsDecomposition(c, x * c, sqrt(b), b / sum(b), cumsum(b) / sum(b))
end

function iterative_pca{T <: Real}(x::Matrix{T})
  error("not yet implemented")
end

function pca{T <: Real}(x::Matrix{T}, method::Symbol)
  if method == :eigenvectors
    eigenvector_pca(x)
  else
    if method == :svd
      svd_pca(x)
    else
      iterative_pca(x)
    end
  end
end

function pca{T <: Real}(x::Matrix{T})
  pca(x, :eigenvectors)
end
