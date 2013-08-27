module DimensionalityReduction
    import Base.print, Base.show, Base.repl_show
    export pca, pcaeig, pcasvd, ica, nmf, mds

    include("types.jl")
    include("pca.jl")
    include("nmf.jl")
    include("mds.jl")
    include("ica.jl")
end
