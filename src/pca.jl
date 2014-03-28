
# Subtracts each column's mean (if center=true),
# divides by each column's standard deviation (if scale=true).
# Returns (scaledData, mean, std), where mean or std may be 
#  empty matrices if center or scale are false
function normalize{T}(X::Matrix{T} ; center=true, scale=true)
    n = size(X,1)

    local m
    if center
        m = mean(X,1)
        X = X .- m
    else
        m = similar(X, 0, 0)
    end

    local s
    if scale
        s = std(X,1)
        X = X ./ s
    else
        s = similar(X, 0, 0)
    end

    return X, m, s
end

# normalizes each column of X, in-place
# returns X itself
function make_cols_unit_norm!( X::Matrix )
    for i=1:size(X,2)
        X[:,i] /= norm(X[:,i])
    end

    return X
end

function pcaeig(X::Matrix ; center=true, scale=true)
    n = size(X,1)
    X, Xmean, Xstd = normalize(X ; center=center, scale=scale)

    C = conj(X' * X / (n - 1))
    XtXeig = eigfact!( C )

    # Eigenvectors are in reverse order from R's
    L = reverse(XtXeig[:values])
    # Zero eigenvalues could have a negative sign
    L = clamp(L, 0.0, Inf)*(n-1)/n
    V = fliplr(XtXeig[:vectors])

    # re-scale back principal components if scaling used
    rotation = scale ? make_cols_unit_norm!(V .* Xstd') : V
    return PCA(rotation, X*V, sqrt(L), L / sum(L), cumsum(L) / sum(L))
end

function pcasvd(X::Matrix ; center=true, scale=true)
    n = size(X,1)
    X, Xmean, Xstd = normalize(X ; center=center, scale=scale)

    Xsvd = svdfact(X)
    pcsd = Xsvd[:S]/sqrt(n)
    pcv = pcsd.^2
    pcvsum = sum(pcv)

    # re-scale back principal components if scaling used
    rotation = scale ? make_cols_unit_norm!(Xsvd[:V] .* Xstd') : Xsvd[:V]
    return PCA(rotation, Xsvd[:U]*Diagonal(Xsvd[:S]), pcsd, pcv/pcvsum, cumsum(pcv)/pcvsum)
end


# TODO: PCA of a DataMatrix, DataArray, DataFrame
#       Should be robust to NA
#       Biplot
#       Eigenvectors orientation? (pcaeig)
#       Summary


# Default uses SVD decomposition
pca(X::Matrix ; kwargs...) = pcasvd(X ; kwargs...)

pcaiterative(X::Matrix) = error("not yet implemented")


