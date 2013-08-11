type PCA
    rotation::Matrix{Float64}
    scores::Matrix{Float64}
    standard_deviations::Vector{Float64}
    proportion_of_variance::Vector{Float64}
    cumulative_variance::Vector{Float64}
end

type ICA
    S::Matrix{Float64}
    H::Matrix{Float64}
    iter::Vector{Int}
    W::Matrix{Float64}
end

type NMF
    W::Matrix{Float64}
    H::Matrix{Float64}
    iterations::Int
    accuracy::Float64
end


function show(io::IO, pc::PCA)
    println(io, "Rotation:")
    show(io, pc.rotation)
    print(io, "\n\n")
    println(io, "Scores:")
    show(io, pc.scores)
    print(io, "\n\n")
    println(io, "Reconstruction:")
    show(io, pc.scores * pc.rotation')
    print(io, "\n\n")
    println(io, "Standard Deviations:")
    show(io, pc.standard_deviations)
    print(io, "\n\n")
    println(io, "Proportion of Variance:")
    show(io, pc.proportion_of_variance)
    print(io, "\n\n")
    println(io, "Cumulative Variance:")
    show(io, pc.cumulative_variance)
end

function show(io::IO, ic::ICA)
    println(io, "Source Matrix:")
    show(io, ic.S)
    print(io, "\n\n")
    println(io, "Mixing Matrix:")
    show(io, ic.H)
    print(io, "\n\n")
    println(io, "Iterations:")
    show(io, ic.iter)
    print(io, "\n\n")
    println(io, "Extracting Matrix:")
    show(io, ic.H)
end


function show(io::IO, res::NMF)
    println(io, "W:")
    show(io, res.W)
    print(io, "\n\n")
    println(io, "H:")
    show(io, res.H)
    print(io, "\n\n")
    println(io, "Reconstruction:")
    show(io, res.W * res.H)
    print(io, "\n\n")
    println(io, "Iterations:")
    show(io, res.iterations)
    print(io, "\n\n")
    println(io, "Accuracy:")
    show(io, res.accuracy)
end
