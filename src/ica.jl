# Translated from Vicente Zarzoso's Matlab code with his permission to use an
# MIT license. Vicente asked us to make clear that this implementation
# has not been reviewed by him and that he cannot be held liable for any
# (a) discrepancies between his code and this code or (b) any consequences
# resulting from any discrepancies between his code and this code.

function ica(X::Matrix; kurtsign=[], tol=.001, max_it=1000, prewhi=false, deftype='r', dimred=false, verbose=false)
    # Kurtosis-based RobustICA method for deflationary ICA/BSS (see references below for details).
    #
    #
    # OUTPUTS:
    #         S       : estimated sources signals (one row per signal, one column per sample)
    #
    #         H       : estimated mixing matrix
    #
    #         iter    : number of iterations (one element per extracted source)
    #
    #         W       : estimated extracting vectors
    #                  (acting on whitened observations if prewhitened is required;
    #                   otherwise, acting on given observations).
    #
    #
    # INPUTS:
    #         X       : observed signals (one row per signal, one column per sample)
    #
    #         kurtsign: source kurtosis signs (one element per source or empty);
    #                   maximize absolute normalized kurtosis if element = 0;
    #                   if empty, maximize absolute normalized kurtosis for all sources
    #
    #         tol     : threshold for statistically-significant termination test of the type
    #                         ||wn - p*w||/||w|| < tol/sqrt(sample size);   (up to a phase shift p)
    #                   termination is also tested by comparing the gradient norm according to:
    #                         ||g|| < tol/sqrt(sample size);
    #                   termination test is not used if tol < 0, so that the algorithm runs the maximum
    #                   number of iterations (except if optimal step size is null)
    #
    #         max_it  : maximum number of iterations per extracted source
    #
    #         prewhi  : prewhitening (via SVD of the observed data matrix);
    #
    #         deftype : deflation type: 'o'rthogonalization, 'r'egression
    #
    #         dimred  : dimensionality reduction in regression if parameter different from zero;
    #                   (not used in deflationary orthogonalization)
    #
    #         verbose : default: quiet operation.
    #
    #
    # REFERENCES:
    #
    # - V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust independent component analysis by iterative maximization</a>
    #   <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the kurtosis contrast with algebraic optimal step size"</a>,
    #   IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
    #
    # - V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/ica07.pdf">"Comparative Speed Analysis of FastICA"</a>,
    #   in: Proceedings ICA-2007, 7th International Conference on Independent Component Analysis
    #       and Signal Separation, London, UK, September 9-12, 2007, pp. 293-300.
    #
    # - V. Zarzoso, P. Comon and M. Kallel,  <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco06.pdf">"How Fast is FastICA?"</a>,
    #   in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference,
    #       Florence, Italy, September 4-8, 2006.
    #
    #
    # HISTORY:
    #   2013-08-01: Conversion from Vicente Zarzoso's Matlab functions by Chris Said.

    Base.depwarn("ica is deprecated.", :ica)

    n, T = size(X)
    X = X - mean(X,2)*ones(1,T)

    # Prewhitening
    if prewhi
        U, D, V = svd(X)
        B = U*D/sqrt(T)     # PCA mixing-matrix estimate
        Z = sqrt(T)*V'      # PCA source estimate
    else
        Z = X;
    end


    # RobustICA algorithm
    dimobs = n;     # number of remaining observations (may change under dimred)
    W = zeros(n,n); # extracting vectors
    I = eye(n);
    P = I;          # proj mtrx for deflationary orthogonalization (if required)
    S = zeros(n,T)  # initializing source matrix

    tol = tol/sqrt(T);         # a statistically-significant termination threshold
    tol2 = sign(tol)*tol^2/2;  # the same thresh in terms of extracting vectors' abs scalar product
    iter = convert(Array{Int}, zeros(n));        # number of iterations

    if isempty(kurtsign)
        kurtsign = zeros(1, n)
    end

    if deftype == 'r'
        do_reg = true;
    elseif deftype == 'o'
        do_reg = false;
        if dimred
            print("Warning: Dimensionality reduction option will be ignored because deflationary orthogonalization has been requested.")
        end
    else
        error("Deflation method not recognized ('deftype' parameter):
                use 'r' or 'o'. Please type 'help(ica)' for details.")
    end



    # iterate over all sources
    for k = 1:n

        it = 0;
        keep_going = true;

        w = I[:, k];    # initialization

        if do_reg
            # keep only required number of components
            w = w[(n-dimobs+1):n]
        end;

        w = w/norm(w);  # normalization

        # project onto extracted vectors' orthog subspace (if defl orthogonalization)
        w = P*w;

        # kurtosis sign of next source to be estimated
        this_kurtsign = kurtsign[k];

        # iterate to extract one source
        while keep_going

            it = it + 1;

            # compute KM optimal step size for gradient descent
            g, mu_opt, norm_g = kurt_gradient_optstep(w, Z, this_kurtsign, P)

            # update extracting vector and project if required
            wn = P*(w + mu_opt*g)

            # normalize
            wn = wn/norm(wn)

            # extracting vector convergence test
            th = abs(1 - abs(wn'*w))[1]

            w = wn;

            if th < tol2 || norm_g < tol || it >= max_it || mu_opt == 0
                # finish when extracting vector converges, the gradient is too small,
                # too many iterations have been run, or the optimal step-size is zero
                keep_going = false;
            end # if th < tol2 ...

        end # while keep_going

        if do_reg
            # zero-padding to account for dimensionality reduction in regression
            W[:,k] = [w; zeros(n-dimobs, 1)]
        else
            # estimated extracting vector
            W[:, k] = w
        end # if do_reg

        s = w'*Z;       # estimated source
        S[k, :] = s;
        iter[k] = it;   # number of  iterations

        if do_reg
            # regression + subtraction
            Z = deflation_regression(Z, s, dimred);

            # recompute observation dimension, in case it's changed during deflation
            dimobs = size(Z, 1)

            # P is not required, but its dimensions should decrease according to Z
            P = eye(dimobs)
        else
            # projection matrix for orthogonalization (if required)
            P = I - W*W'
        end # if do_reg

    end # for k

    # Mixing matrix estimation. LS estimate.
    H = X*S'*pinv(S*S');


    return ICA(S, H, iter, W)
end


function kurt_gradient_optstep(w, X, s, P)

    # Computes optimal step size in the gradient-based optimization of the normalized kurtosis contrast
    # (single iteration).
    #
    #
    # OUTPUTS:
    #         g      : search direction (normalized gradient vector)
    #
    #         mu_opt : optimal step size globally optimizing the normalized kurtosis contrast function
    #                  along direction g from f
    #
    #         norm_g : non-normalized gradient vector norm.
    #
    # INPUTS:
    #         w      : current extracting vector coefficients
    #
    #         X      : sensor-output data matrix (one signal per row, one sample per column)
    #
    #         s      : source kurtosis sign; if zero, the maximum absolute value of the contrast is sought
    #
    #         P      : projection matrix (used in deflationary orthogonalization; identity matrix otherwise).



    L, T = size(X);

    mu_opt = 0;     # default optimal step-size value
    norm_g = 0;     # initialize gradient norm


    # Compute search direction (gradient vector)

    y = w'*X

    ya2 = y.*conj(y);
    y2 = y.*y;
    ya4 = ya2.*ya2;

    Eya2 = mean(ya2);
    Ey2 = mean(y2);
    Eya4 = mean(ya4);

    #check for zero denominator
    if abs(Eya2) < eps()
        g = zeros(L, 1);
        norm_g = 0;
    else

        # compute gradient if contrast denominator is not null
        Eycx = X*y'/T
        Eyx = X*y.'/T
        Ey3x = X*(ya2.*y)'/T

        # contrast numerator and denominator at current point
        p1 = Eya4 - abs(Ey2)^2
        p2 = Eya2

        g = 4*( (Ey3x - Eyx*Ey2')*p2 - p1*Eycx )/p2^3

        # project if required (normalize later)
        g = P*g

        norm_g = norm(g)

        if norm_g >= 10e-13 # eps() in original Matlab implementation

            # normalize the gradient
            g = g/norm_g

            # Compute optimal step size
            gg = g'*X


            # calculate interim values for contrast rational function
            ya2 = y.*conj(y);
            ga2 = gg.*conj(gg);
            reygc = real(y.*conj(gg));

            g2 = gg.*gg;
            yg = y.*gg;

            Eya2reygc = mean(ya2.*reygc)
            Ereygc2 = mean(reygc.^2)
            Ega2reygc = mean(ga2.*reygc)
            Ega4 = mean(ga2.^2)
            Eya2ga2 = mean(ya2.*ga2)

            Ega2 = mean(ga2)
            Ereygc = mean(reygc)

            Eg2 = mean(g2)
            Eyg = mean(yg)

            h0 = Eya4 - abs(Ey2)^2
            h1 = 4*Eya2reygc - 4*real(Ey2*Eyg')
            h2 = 4*Ereygc2 + 2*Eya2ga2 - 4*abs(Eyg)^2 - 2*real(Ey2*Eg2')
            h3 = 4*Ega2reygc - 4*real(Eg2*Eyg')
            h4 = Ega4 - abs(Eg2)^2

            P = [h4 h3 h2 h1 h0]

            i0 = Eya2
            i1 = 2*Ereygc
            i2 = Ega2

            Q = [i2 i1 i0]

            # normalized kurtosis contrast = P/Q^2 - 2

            a0 = -2*h0*i1 + h1*i0;
            a1 = -4*h0*i2 - h1*i1 + 2*h2*i0;
            a2 = -3*h1*i2 + 3*h3*i0;
            a3 = -2*h2*i2 + h3*i1 + 4*h4*i0;
            a4 = -h3*i2 + 2*h4*i1;

            p = [a4 a3 a2 a1 a0]



            #----- Find the roots of p -----#

            n = length(p)

            inz = find(p) #indices not zero


            # Strip leading zeros and throw away.
            # Strip trailing zeros, but remember them as roots at zero.
            nnz = length(inz)

            p = p[inz[1]:inz[nnz]]
            r = zeros(n-inz[nnz],1)

            # Prevent relatively small leading coefficients from introducing Inf
            # by removing them.
            d = p[2:end]./p[1]
            while any(isinf(d))
                c = p[2:end]
                d = p[2:end]./p[1]
            end

            # Polynomial roots via a companion matrix
            n = length(p)
            if n > 1
                a = zeros(n-1,n-1)
                for i = 1:(n-2)
                    a[i+1,i] = 1
                end
                a[1,:] = -d
                evals, evecs = eig(a)
                r = [r; evals]
            end

            #-------------------------------#


            rr = real(r);             # keep real parts only
            Pval = zeros(length(rr))  # initializing
            Q2val = zeros(length(rr)) # initializing

            for i = 1:length(rr)
                Pval[i] = polyval(P, rr[i])
                Q2val[i] = polyval(Q, rr[i])^2
            end


            # check roots not shared by denominator
            # NOTE: in theory, the denominator can never
            # cancel out if the gradient is used as a search direction,
            # due to the orthogonality between the extracting vector and the
            # corresponding gradient. (Only exception: if it is the last source
            # to be extracted. But this scenario is detected by the gradient norm)
            nonzero_Q2val = find(Q2val .> eps())


            if !isempty(nonzero_Q2val)
                Pval = Pval[nonzero_Q2val];
                Q2val = Q2val[nonzero_Q2val];
                rr = rr[nonzero_Q2val];

                Jkm_val = Pval./Q2val - 2;    # normalized kurtosis


                if s != 0
                    # maximize or minimize kurtosis value, depending on kurt sign
                    Jkm_val = real(s*Jkm_val);
                else
                    # maximize absolute kurtosis value, if no sign is given
                    Jkm_val = abs(Jkm_val);
                end

                Jmax, im = findmax(Jkm_val);
                mu_opt = rr[im];              # optimal step size

            end # if isempty(nonzero_Q2val)

        end # if norm(g) < eps()

    end # if abs(Eya2) < eps()

    return g, mu_opt, norm_g
end

function polyval(p, x)
    n = length(p) - 1
    y = 0
    for i = 0:n
        y = y + p[n+1-i]*x^i
    end
    return y
end



function deflation_regression(Z, s, dimred);

    # Performs deflation by subtracting the estimated source contribution to the observations as:
    #
    #    Z' = Z - h*s
    #
    # The source direction h is estimated via the least squares solution to the linear regression
    # problem:
    #
    #    h_opt = arg min_h ||Z - h*s||^2 = Z*s'/(s*s').
    #
    #
    #
    # OUTPUT:
    #         Zn     : observed data after subtraction (one signal per row, one sample per column).
    #
    #
    # INPUTS:
    #         Z      : observed data (one signal per row, one sample per column)
    #
    #         s      : estimated source (row vector with one sample per column)
    #
    #         dimred : perform dimensionality reduction if parameter different from zero.



    # extracted source power times sample size
    s2 = s*s';


    if abs(s2)[1] > eps()

        # source direction estimated via least squares
        h = Z*s'/s2;

        if dimred

            n = length(h)
            Q = [h eye(n, n-1)]
            Q, R = qr(Q)

            # orthonormal basis of orthogonal subspace of h
            Q = Q[:,2:n]

            # remaining contribution with dimensionality reduction
            Z = Q'*Z;

        else

            # without dimensionality reduction
            Z = Z - h*s;

        end # if dimred

    end # if s2 > eps

    return Z
end
