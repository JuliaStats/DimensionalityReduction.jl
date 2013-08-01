using DimensionalityReduction


# Function compares estimated sources to true sources
function evaluate_source_estimation(S, S_true)

    # The order of estimated sources will differ from the arbitrary order
    # of true sources. Here we permute estimated sources to the order of
    # the true sources using a method that should work fine for well-separated
    # sources.
    P = zeros(K, K)        # permutation matrix
    for k = 1:K
        s = S[k, :]
        r = S_true*s'      # cross-correlation
        pos_max = findmax(abs(r))[2]
        P[pos_max, k] = (r[pos_max]/(s*s'))[1]  # ordering with optimum scale
    end

    S = P*S

    SMSE = sum(sum(abs((S_true - S).^2)))/(K*T)

    return SMSE
end


T = 1000          # sample size
K = 10            # number of sources and observations

# true source matrix
S_true = zeros(K, T);

# a few sub-Gaussian real-valued sources
sub_g = (rand(5, T) .> 0.5)

# a few super-Gaussian real-valued sources
sup_g = (rand(5, T) .> 0.9)

# create source matrix
S_true[[sub_g; sup_g]] = 1

# Normalize the sources
S_true = S_true .- mean(S_true,2)    # remove mean
S_true = S_true ./ std(S_true,2)     # unit-variance normalization


# generate random mixing matrix
H_true = randn(K, K)

# generate observed signals
X = H_true*S_true

# Evaluate ICA algorithm with different parameter cominations

results = ica(X, kurtsign=[], tol=.001, max_it=1000, prewhi=true, deftype='r', dimred=false, verbose=false)
SMSE = evaluate_source_estimation(results.S, S_true)
@assert SMSE < .1


results = ica(X, kurtsign=[], tol=.001, max_it=1000, prewhi=true, deftype='o', dimred=false, verbose=false)
SMSE = evaluate_source_estimation(results.S, S_true)
@assert SMSE < .1


results = ica(X, kurtsign=[], tol=.001, max_it=1000, prewhi=false, deftype='r', dimred=false, verbose=false)
SMSE = evaluate_source_estimation(results.S, S_true)
@assert SMSE < .1


#Lots of iterations (>> 10000) are required for this parameter combination to perform well.
results = ica(X, kurtsign=[], tol=.001, max_it=1000, prewhi=false, deftype='o', dimred=false, verbose=false)
SMSE = evaluate_source_estimation(results.S, S_true)
@assert SMSE < 4


results = ica(X, kurtsign=[], tol=.001, max_it=1000, prewhi=true, deftype='r', dimred=true, verbose=false)
SMSE = evaluate_source_estimation(results.S, S_true)
@assert SMSE < .1


results = ica(X, kurtsign=[], tol=.001, max_it=1000, prewhi=false, deftype='r', dimred=true, verbose=false)
SMSE = evaluate_source_estimation(results.S, S_true)
@assert SMSE < .1
