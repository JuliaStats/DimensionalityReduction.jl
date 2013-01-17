function pcaeig(X::Matrix)
    L, Z = eig(cov(X))
    P = Z'
    # Eigenvectors are in reverse order from R's
    # Need to clamp any negative eigenvalues before calling sqrt()
    L = reverse(L)
    for i in 1:length(L)
        L[i] = clamp(L[i], 0.0, Inf)
    end
    P = fliplr(P)
    Y = X * P # X ~ Y * P'
    return PCA(P, Y, sqrt(L), L / sum(L), cumsum(L) / sum(L))
end

function pcasvd(X::Matrix)
    A, B, C = svd(cov(X))
    return PCA(C, X * C, sqrt(B), B / sum(B), cumsum(B) / sum(B))
end

pcaiterative(x::Matrix) = error("not yet implemented")

# TODO: PCA of a DataMatrix
#       Should be robust to NA

pca(x::Matrix) = pcasvd(x)
