function pcaeig(X::Matrix)
    n = size(X,1)
    XtXeig = eigfact!(cov(X))
    # Eigenvectors are in reverse order from R's
    # Need to clamp any negative eigenvalues before calling sqrt()
    L = reverse(XtXeig[:values])
    for i in 1:length(L)
        L[i] = clamp(L[i], 0.0, Inf)*(n-1)/n
    end
    Z = fliplr(XtXeig[:vectors])
    return PCA(Z, X*Z, sqrt(L), L / sum(L), cumsum(L) / sum(L))
end

function pcasvd(X::Matrix)
    Xsvd = svdfact(X)
    pcsd = Xsvd[:S]/sqrt(size(X,1))
    pcv = pcsd.^2
    pcvsum = sum(pcv)
    return PCA(Xsvd[:V], Xsvd[:U]*Diagonal(pcsd), pcsd, pcv/pcvsum, cumsum(pcv)/pcvsum)
end

pcaiterative(x::Matrix) = error("not yet implemented")

# TODO: PCA of a DataMatrix
#       Should be robust to NA

pca(x::Matrix) = pcasvd(x)
