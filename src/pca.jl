
# Subtracts each column's mean (if center=true),
# divides by each column's standard deviation (if scale=true).
function normalize(X::Matrix ; center=true, scale=true)
    n = size(X,1)
    if center
        X = X.-mean(X,1)
    end
    if scale
        X = X./std(X,1)
    end
    return X
end

function pcaeig(X::Matrix ; center=true, scale=true)
    n = size(X,1)
    X = normalize(X ; center=center)
    C = scale ? cor(X) : cov(X)
    XtXeig = eigfact!(C)
    # Eigenvectors are in reverse order from R's
    L = reverse(XtXeig[:values])
    # Zero eigenvalues could have a negative sign
    L = clamp(L, 0.0, Inf)*(n-1)/n
    V = fliplr(XtXeig[:vectors])
    return PCA(V, X*V, sqrt(L), L / sum(L), cumsum(L) / sum(L))
end

function pcasvd(X::Matrix ; center=true, scale=true)
    n = size(X,1)
    X = normalize(X ; center=center, scale=scale)
    Xsvd = svdfact(X)
    pcsd = Xsvd[:S]/sqrt(n)
    pcv = pcsd.^2
    pcvsum = sum(pcv)
    return PCA(Xsvd[:V], Xsvd[:U]*Diagonal(Xsvd[:S]), pcsd, pcv/pcvsum, cumsum(pcv)/pcvsum)
end


# TODO: PCA of a DataMatrix, DataArray, DataFrame
#       Should be robust to NA
#       Biplot
#       Eigenvectors orientation? (pcaeig)
#       Summary


# Default uses SVD decomposition
pca(X::Matrix ; kwargs...) = pcasvd(X ; kwargs...)

pcaiterative(X::Matrix) = error("not yet implemented")


