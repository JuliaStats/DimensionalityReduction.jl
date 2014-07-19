module DimensionalityReduction
	warn("The DimensionalityReduction package is deprecated. It is superseded by a new package MultivariateStats.")

    export pca, pcaeig, pcasvd
    export ica, nmf, mds, tsne

    include("types.jl")
    include("pca.jl")
    include("nmf.jl")
    include("mds.jl")
    include("ica.jl")
    include("tsne.jl")
end
