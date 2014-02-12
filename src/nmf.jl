# Reference
#
# "Algorithms for Non-negative Matrix Factorization"
# by Daniel D Lee, Sebastian H Seung

# Numeric stability ideas borrowed from the Milk implementation of NMF
# Convergence diagnostic borrowed from Milk, although it seems unhelpful
# TODO: Reimplement my randomized NMF and compare with this version

function nmf(X::Matrix,
             d::Integer,
             max_iter::Integer,
             tolerance::FloatingPoint,
             epsilon::FloatingPoint)

    Base.depwarn("nmf is deprecated. Please checkout the NMF package.", :nmf)

    n, m = size(X)

    W = randn(n, d).^2
    H = randn(d, m).^2

    iter = 0
    while iter < max_iter
        iter += 1
        updateH = (W' * X) ./ (W' * W * H + epsilon)
        H .*= updateH
        updateW = (X * H') ./ (W * (H * H') + epsilon)
        W .*= updateW
        max_update = max(maximum(updateH), maximum(updateW))
        if abs(1.0 - max_update) < tolerance
            break
        end
    end

    return NMF(W, H, iter, norm(X - W * H))
end

nmf(X::Matrix, d::Integer) = nmf(X, d, 1_000, 10e-8, 10e-8)
