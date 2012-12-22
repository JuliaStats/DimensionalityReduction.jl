require("DataFrames")
using DataFrames

module DimensionalityReduction
	using DataFrames

	import Base.print, Base.show, Base.repl_show
	export pca, ica, nmf

	require("DimensionalityReduction/src/types.jl")
	require("DimensionalityReduction/src/pca.jl")
end
