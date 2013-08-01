module DimensionalityReduction
	import Base.print, Base.show, Base.repl_show
	export pca, pcaeig, pcasvd, ica, nmf, mds

	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "types.jl"))
	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "pca.jl"))
	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "nmf.jl"))
	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "mds.jl"))
    include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "ica.jl"))
end
