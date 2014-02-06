module DimensionalityReduction
    export pca, pcaeig, pcasvd, biplot, pca_example
    export ica, nmf, mds, tsne

    include("types.jl")
    include("pca.jl")
    include("nmf.jl")
    include("mds.jl")
    include("ica.jl")
    include("tsne.jl")
end
