module DimensionalityReduction
    using Distance

    export pca, pcaeig, pcasvd
    export ica, nmf, mds, sne, tsne

    include("types.jl")
    include("pca.jl")
    include("nmf.jl")
    include("mds.jl")
    include("ica.jl")
    include("SNE/PQ.jl")
    include("SNE/cost.jl")
    include("SNE/fit.jl")
end
