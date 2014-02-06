

# Subtracts each column's mean (if center=true),
# divides by each column's standard deviation (if scale=true).
function normalize!(X::Matrix ; center=true, scale=false)
    n = size(X,1)
    if center
        means = mean(X,1)
        for i in 1:n
            X[i,:] = X[i,:].-means
        end
    end
    if scale
        sd = std(X,1)
        for i in 1:n
            X[i,:] = X[i,:]./sd
        end
    end
    return X
end
normalize(X::Matrix ; kwargs...) = normalize!(copy(X) ; kwargs...)

#-----------------------------------------------------------------------------#

function pcaeig(X::Matrix ; correlation=false)
    n = size(X,1)
    C = correlation ? cor(X) : cov(X)
    XtXeig = eigfact!(C)
    # Eigenvectors are in reverse order from R's
    L = reverse(XtXeig[:values])
    V = fliplr(XtXeig[:vectors])
    return PCA(V, X*V, sqrt(L), L / sum(L), cumsum(L) / sum(L))
end

function pcasvd(X::Matrix ; center=true, scale=false, tol=Inf)
    normalize!(X ; center=center, scale=scale)
    n = size(X,1)
    Xsvd = svdfact(X)
    pcsd = Xsvd[:S]/sqrt(n)
    pcv = pcsd.^2
    pcvsum = sum(pcv)
    return PCA(Xsvd[:V], Xsvd[:U]*Diagonal(Xsvd[:S]), pcsd, pcv/pcvsum, cumsum(pcv)/pcvsum)
end

#-----------------------------------------------------------------------------#

# Not functional (where to do the import of Gadfly??)
function biplot(P::PCA)
    using Gadfly
    scores = P.scores[:,1:2]
    return plot(x=scores[1],y=scores[2], Geom.point)
end

# Not functional (see above)
function pca_example()
    using RDatasets
    df = data("datasets", "iris")
    iris = convert(Array,DataArray(df[:,1:4]))
    psvd = pca(iris)
    pl = biplot(psvd)
end

#-----------------------------------------------------------------------------#

# TODO: PCA of a DataMatrix, DataArray, DataFrame
#       Should be robust to NA
#       Biplot
#       Eigenvectors orientation

#-----------------------------------------------------------------------------#

# Default uses SVD decomposition
pca(x::Matrix ; kwargs...) = pcasvd(x ; kwargs...)

pcaiterative(x::Matrix) = error("not yet implemented")


