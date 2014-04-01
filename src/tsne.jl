function Hbeta!(P::Vector,
                D::Array,
                beta::Real = 1.0)
    # Compute the perplexity and the P-row for a specific value of 
    # the precision of a Gaussian distribution.
    n = length(D)
    for i in 1:n
        P[i] = exp(-D[i] * beta)
    end
    sumP, sumDP = 0.0, 0.0
    for i in 1:n
        sumP += P[i]
        sumDP += D[i] * P[i]
    end
    for i in 1:length(P)
        P[i] = P[i] / sumP
    end
    H = log(sumP) + beta * sumDP / sumP
    return H
end
    
function x2p(X::Matrix,
             tol::Real = 1e-5,
             perplexity::Real = 30.0,
             tracing::Bool = false)
    # Performs a binary search to get P-values in such a way that
    # each conditional Gaussian has the same perplexity.
    n, d = size(X)
    sum_X = sum(X.^2, 2)
    D = (-2 * (X * X') .+ sum_X)' .+ sum_X 
    P = zeros(n, n)
    beta = ones(n, 1)
    logU = log(perplexity)

    # Loop over all datapoints
    range = [1:n]
    thisP = Array(Float64, n - 1)
    for i in 1:n
        # Print progress
        if tracing
            if mod(i, 500) == 0
                @printf "Computing P-values for point %d of %d...\n" i n
            end
        end

        # Compute the Gaussian kernel and entropy for the current precision
        betamin = -Inf
        betamax = Inf

        inds = range[range .!=i]
        Di = D[i, inds]
        H = Hbeta!(thisP, Di, beta[i])

        # Evaluate whether the perplexity is within tolerance
        Hdiff = H - logU
        tries = 0
        while abs(Hdiff) > tol && tries < 50
            # If not, increase or decrease precision
            if Hdiff > 0
                betamin = beta[i]
                if isinf(betamax)
                    beta[i] *= 2
                else
                    beta[i] = (beta[i] + betamax) / 2
                end
            else
                betamax = beta[i];
                if isinf(betamin)
                    beta[i] /= 2
                else
                    beta[i] = (beta[i] + betamin) / 2
                end
            end

            # Recompute the values
            H = Hbeta!(thisP, Di, beta[i])
            Hdiff = H - logU
            tries += 1
        end
        # Set the final row of P
        P[i, inds] = thisP
    end

    if tracing
        @printf "Mean value of sigma: %f\n" mean(sqrt(1 / beta))
    end

    # Return final P-matrix
    return P
end

function innerpca{T}(X::Matrix{T}, ndims::Integer = 50)
    # Runs PCA on the NxD array X in order to reduce its
    # dimensionality to ndims dimensions.
    n, d = size(X)
    X = X .- mean(X, 1)
    Xp = float(X' * X)
    M = eigfact(Xp).vectors::typeof(Xp)
    Y = X * M[:, 1:min(d, ndims)]
    return Y
end

function tsne(X::Matrix,
              no_dims::Integer = 2,
              initial_dims::Integer = -1,
              max_iter::Integer = 1000,
              perplexity::Real = 30.0,
              eta::Real = 500.0,
              min_gain::Real = 0.01,
              initial_momentum::Real = 0.5,
              final_momentum::Real = 0.8,
              tracing::Bool = false)
    #Runs t-SNE on the dataset in the NxD array X to reduce its dimensionality to no_dims dimensions.

    Base.depwarn("tsne is deprecated.", :tsne)

    # Initialize variables
    if initial_dims > 0
        X = innerpca(X, initial_dims)
    end

    n, d = size(X)

    Y = randn(n, no_dims)
    dY = zeros(n, no_dims)
    iY = zeros(n, no_dims)
    gains = ones(n, no_dims)
    
    # Compute P-values
    P = x2p(X, 1e-5, perplexity, tracing)
    nP = n * n
    # P = P + P'
    for j in 1:n
        for i in j:n
            pij, pji = P[i, j], P[j, i]
            pijnew = pij + pji
            P[i, j] = pijnew
            P[j, i] = pijnew
        end
    end
    # P = P / sum(P)
    # P = P * 4
    sumP = sum(P)
    for i in 1:nP
        P[i] /= sumP
        P[i] *= 4 # early exaggeration
    end
    # NO-OP
    # P = maximum(P, 1e-12)

    # Run iterations
    sum_Y = Array(Float64, n)
    Q = Array(Float64, n, n)
    PQ = Array(Float64, n, n)
    logs = Array(Float64, n, n)
    YYt = Array(Float64, n, n)
    inter = Array(Float64, n, n)
    num = Array(Float64, n, n)
    tmp = Array(Float64, n, no_dims)

    for iter in 1:max_iter
        # Compute pairwise affinities as sum of squared row entries
        for i in 1:n
            sum_Y[i] = 0.0
        end
        for j in 1:no_dims
            for i in 1:n
                sum_Y[i] += Y[i, j]^2
            end
        end

        A_mul_Bt!(YYt, Y, Y)

        # Compute ((-2 * YYt) .+ sum_Y)
        for j in 1:n
            for i in 1:n
                inter[i, j] = -2 * YYt[i, j] + sum_Y[i]
            end
        end

        # Compute 1 / (1 + ((-2 * (YYt)) .+ sum_Y)' .+ sum_Y)
        for j in 1:n
            for i in 1:n
                num[i, j] = 1 / (1 + inter[j, i] + sum_Y[i])
            end
        end

        # Set diagonal to zero
        for i in 1:n
            num[i, i] = 0.0
        end

        sum_num = sum(num)
        for i in 1:length(Q)
            Q[i] = num[i] / sum_num
        end
        # NO-OP
        # Q = maximum(Q, 1e-12)

        # Compute gradient
        for i in 1:length(PQ)
            PQ[i] = P[i] - Q[i]
        end
 
        for k in 1:n
            for j in 1:no_dims
                for i in 1:n
                    tmp[i, j] = PQ[i, k] * num[i, k]
                end
            end
            for j in 1:no_dims
                for i in 1:n
                    tmp[i, j] = tmp[i, j] * (Y[k, j] - Y[i, j])
                end
            end
            for j in 1:no_dims
                s = 0.0
                for i in 1:n
                    s += tmp[i, j]
                end
                dY[k, j] = s
            end
        end

        # Perform the update
        momentum = iter <= 20 ? initial_momentum : final_momentum
        for i in 1:length(gains)
            gains[i] = (gains[i] + 0.2) * ((dY[i] > 0) != (iY[i] > 0)) +
                       (gains[i] * 0.8) * ((dY[i] > 0) == (iY[i] > 0))
            if gains[i] < min_gain
                gains[i] = min_gain
            end
        end
        for i in 1:length(iY)
            iY[i] = momentum * iY[i] - eta * (gains[i] * dY[i])
            Y[i] = Y[i] + iY[i]
        end

        # Subtract column means from Y
        for j in 1:no_dims
            m = 0.0
            for i in 1:n
                m += Y[i, j]
            end
            m /= n
            for i in 1:n
                Y[i, j] -= m
            end
        end

        # Compute current value of cost function
        if mod(iter + 1, 10) == 0
            for i in 1:length(logs)
                logs[i] = log(P[i] / Q[i])
                # So we don't get NaN when the error is computed
                logs[i] = isnan(logs[i]) ? 0.0 : logs[i]
            end
            C = 0.0
            for i in 1:nP
                C += P[i] * logs[i]
            end
            if tracing
                @printf "Iteration %d: error %f\n" iter + 1 C
            end
        end

        # Restore P to its proper magnitude
        if iter == 100
            for i in 1:nP
                P[i] = P[i] / 4
            end
        end
    end

    # Return solution
    return TSNE(Y)
end
