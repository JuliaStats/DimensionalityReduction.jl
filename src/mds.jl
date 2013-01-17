function mds(D::Matrix, k::Integer)
    n = size(D, 1)
    P = D.^2
    J = eye(n) - (1.0 / n) * ones(n, n)
    B = -0.5 * J * P * J
    L, Z = eig(B)
    # Only use non-complex, positive eigenvalues
    indices = find(imag(L) != 0.0 && real(L) .> 0)
    L = real(L[indices])
    Z = real(Z[:, indices])
    # If there are no longer k eigenvalues left, shrink k and warn user
    m = length(L)
    if k > m
        warn(@sprintf "%d dimensional MDS not possible. Using %d instead." k m)
        k = m
    end
    indices = (m - k + 1):m
    X = Z[:, indices] * diagm(sqrt(L[indices]))
    return MDS(X, D, k)
end

mds(D::Matrix) = mds(D, 2)
