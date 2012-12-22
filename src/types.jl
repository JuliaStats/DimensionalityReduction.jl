type PrincipalComponentsDecomposition
  rotation::Matrix{Float64}
  scores::Matrix{Float64}
  standard_deviations::Vector{Float64}
  proportion_of_variance::Vector{Float64}
  cumulative_variance::Vector{Float64}
end

function print(io::IOStream, pc::PrincipalComponentsDecomposition)
  println(io, "Principal Components:")
  println(io, pc.rotation)
  println(io)
  println(io, "Proportion of Variance:")
  println(io, pc.proportion_of_variance)
  println(io)
  println(io, "Cumulative Variance:")
  println(io, pc.cumulative_variance)
  println(io)
end
print(pc::PrincipalComponentsDecomposition) = print(OUTPUT_STREAM, pc)

function show(io::IOStream, pc::PrincipalComponentsDecomposition)
  print(io, pc)
end
show(pc::PrincipalComponentsDecomposition) = show(OUTPUT_STREAM, pc)

function repl_show(io::IOStream, pc::PrincipalComponentsDecomposition)
  print(io, pc)
end

repl_show(pc::PrincipalComponentsDecomposition) = repl_show(OUTPUT_STREAM, pc)
