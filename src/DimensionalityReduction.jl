module DimensionalityReduction
	import Base.print, Base.show, Base.repl_show
	export pca, pcaeig, pcasvd, ica, nmf

	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "types.jl"))
	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "pca.jl"))
	include(joinpath(julia_pkgdir(), "DimensionalityReduction", "src", "nmf.jl"))
end
