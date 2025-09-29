# functions to load for ILR method

# R emulation of GAUSS eye(n)
eye <- function(n) {
    return (diag(x = 1, nrow = n))
}

# R emulation of GAUSS maxc(x)
maxc <- function(x) {
    return (matrix(data = apply(as.matrix(x), 2, max), ncol = 1))
}

# R emulation of GAUSS minc(x)
minc <- function(x) {
    return (matrix(data = apply(as.matrix(x), 2, min), ncol = 1))
}

# R emulation of GAUSS sumc(x)
sumc <- function(x) {
    return (matrix(data = apply(as.matrix(x), 2, sum), ncol = 1))
}    

# R emulation of GAUSS vec(x)
vec <- function(x) {
    return (matrix(data = as.matrix(x), nrow = length(x), ncol = 1))
}

# R emulation of GAUSS vecr(x)
vecr <- function(x) {
    return (matrix(data = t(as.matrix(x)), nrow = length(x), ncol = 1))
}

# R emulation of GAUSS reshape(x, r, c)
reshapex <- function(x, m, n) {

    # Get some useful quantities
    x <- as.matrix(x)
    i <- nrow(x)
    j <- ncol(x)

    size_x <- i * j
    size_y <- m * n
    
    if (size_x == size_y) {
        y <- vec(t(x))
        y <- matrix(data = y, nrow = m, ncol = n, byrow = TRUE)
    } else if(size_x < size_y) {
        y <- rep(vec(t(x)), times = size_y)
        y <- matrix(data = y, nrow = m, ncol = n, byrow = TRUE)
    } else {
        y <- vec(t(x))
        y <- matrix(data = y[1:size_y], nrow = m, ncol = n, byrow = TRUE)
    } 
    return (y)
}

# Transform dependent and explanatory variables for estimation
transf <- function(y, z, S) {
    
    # Calculate a few useful quantities
    y <- as.matrix(y)
    z <- as.matrix(z)
    m <- nrow(y)
    n <- ncol(y)
    ID <- eye(n)
    
    maty <- t(y[1, , drop = FALSE])
# browser()
    matz <- ID %x% z[1, , drop = FALSE] %*% S
 
    i <- 2L
    while (i <= m) {
        maty <- rbind(maty, t(y[i, , drop = FALSE]))
        matz <- rbind(matz, ID %x% z[i, , drop = FALSE] %*% S)
        i <- i + 1L
    }
    return (list(maty = maty, matz = matz))
}

# Compute the square root of a symmetric positive definition matrix
sqrm <- function(m) {
    # check dimensions
    if (is.null(dim(m))) {
        sqm <- sqrt(m)
    } else {
        tmp <- eigen(m, symmetric = TRUE)
        D <- diag(x = tmp$values, nrow = nrow(m))
        V <- tmp$vectors
        sqm <- V %*% sqrt(D) %*% t(V)
    }
    return (sqm)
}

# Estimate the coefficients under the optimal partitions
estim <- function(maty, matz, n, m, br, max_iter = 200) {
    bhat <- matrix(data = 0.0, nrow = (m + 1), ncol = ncol(matz))
    vv <- matrix(data = 0.0, nrow = ((m + 1) * n), ncol = n)
    br <- rbind(0, br, (nrow(maty) / n))
    k <- 1L
    while (k <= (m + 1)) {
        i <- br[k] + 1
        j <- br[k + 1]

        segy <- maty[((i - 1)*n+1):(j * n), , drop = FALSE]
        segx <- matz[((i - 1)*n+1):(j * n), , drop = FALSE]
        b <- qr.coef(qr = qr(segx), y = segy)
        res <- segy - segx %*% b
        umat <- reshapex(res, nrow(res)/n, n)
        vvar <- crossprod(x = umat, y = NULL) / (j - i + 1)
        vstar <- vvar + 1
        itr <- 1L
        while ((maxc(abs(vec(t(vvar - vstar)))) > 1E-6) && (itr < max_iter)) {
            vstar <- vvar
            ibigv <- eye(nrow(segy)/n) %x% qr.solve(vstar)
            b <- qr.solve(t(segx) %*% ibigv %*% segx, t(segx) %*% ibigv %*% segy)
            res <- segy - segx %*% b
            umat <- reshapex(res, nrow(res)/n, n)
            vvar <- crossprod(x = umat, y = NULL) / (j - i + 1)
            itr <- itr + 1L
            if (itr == max_iter) {
                warning("FGLS failed to converge\n", call. = FALSE)
                break
            }
        }
        
        bhat[k, ] <- t(b)
        vv[((k - 1) * n + 1):(k * n), ] <- vvar
        k <- k + 1L
    }

    return (list(bhat = bhat, vv = vv))

}

# Estimate the model allowing for cross regime (linear) restrictions
restim <- function(maty, matz, bigt, n, m, br, R, brbeta, brv, max_iter = 1000) {

# Inputs:
#   maty = dependent variable
#   matz = independent variable(s)
#   n = number of equations
#   m = number of breaks
#   bigt = sample size
#   br = break dates
#   R = the restriction, taking the form beta = R %*% teta, where teta is the vector of basic 
#       parameters of the restricted model.
#   brv:
#       brv == 1 if the coefficients are allowed to change
#       brv == 0 otherwise
#   brbeta:
#       brbeta == 1 if the coefficients are allowed to change
#       brbeta == 0 otherwise
#
# Output (list object):
#   nbeta = the estimate of coefficients imposing restrictions
#   nvv = the estimate of covariances imposing restrictions
#
# Note:
#   the model includes partial structural change and block-partial structural change models

seg <- rbind(0, br, (nrow(maty) / n))

# Now estimate an unrestricted model using estim()

tmp <- estim(maty, matz, n, m, br)
beta <- tmp$bhat
vv <- tmp$vv

# Changes in both coefficients and covariances
if ((brv == 1) && (brbeta == 1)) {
    pmatz <- pzbar(matz, m, br, bigt) # %*% R
    pmatz <- pmatz %*% R
    k <- 1L
    ibigv <- eye(bigt * n)
    while (k <= (m + 1)) {
        i <- seg[k] + 1
        j <- seg[(k + 1)]
        ibigv[((i - 1)*n+1):(j * n), ((i - 1)*n+1):(j * n)] <- eye(j - i + 1) %x% qr.solve(vv[((k - 1)*n+1):(k * n), ])
        k <- k + 1L
    }

    b <- qr.solve(t(pmatz) %*% ibigv %*% pmatz, t(pmatz) %*% ibigv %*% maty)
    bstar <- b + 10
    itr <- 1L
    while ((maxc(abs(bstar - b)) > 1E-6) && (itr < max_iter)) {
        bstar <- b
        k <- 1L
        # Update the VCV matrix
        while (k <= (m + 1)) {
            i <- seg[k] + 1
            j <- seg[k + 1]
            tmpy <- maty[((i - 1)*n+1):(j * n), , drop = FALSE]
            tmpx <- pmatz[((i - 1)*n+1):(j * n), , drop = FALSE]
            res <- tmpy - tmpx %*% bstar
            umat <- reshapex(res, nrow(res)/n, n)
            vvar <- crossprod(x = umat, y = NULL) / (j - i + 1)
            vv[((k - 1)*n+1):(k * n), ] <- vvar
            k <- k + 1L
        }
        
        k <- 1L
        # Create the block diagonal matrix of the covariances
        while (k <= (m + 1)) {
            i <- seg[k] + 1
            j <- seg[(k + 1)]
            ibigv[((i - 1)*n+1):(j * n), ((i - 1)*n+1):(j * n)] <- eye(j - i + 1) %x% qr.solve(vv[((k - 1)*n+1):(k * n), ])
            k <- k + 1L
        }
        
        b <- qr.solve(t(pmatz) %*% ibigv %*% pmatz, t(pmatz) %*% ibigv %*% maty)
        itr <- itr + 1L
        if (itr == max_iter) {
            warning("Iteration has reached the upper bound.\n", call. = FALSE)
            break
        }
    }

    nbeta <- R %*% b
    nvv <- vv
}

# Changes in covariances only
# iterative method is used for the optimisation, with the unrestricted estimate as the prior input.
# The procedure is iterated until convergence.
# if ((brv == 1) && (brbeta != 1)) {
#     ibigv <- eye(bigt * n)
#     k <- 1L
#     while (k <= (m + 1)) {
#         i <- seg[k] + 1
#         j <- seg[(k + 1)]
#         ibigv[((i - 1)*n+1):(j * n), ((i - 1)*n+1):(j * n)] <- eye(j - i + 1) %x% qr.solve(vv[((k - 1)*n+1):(k * n), ])
#         k <- k + 1L
#     }

#     b <- qr.solve(t(matz) %*% ibigv %*% matz, t(matz) %*% ibigv %*% maty)
#     bstar <- b + 10
#     itr <- 1L
#     while ((maxc(abs(bstar - b)) > 1E-6) && (itr < max_iter)) {
#         bstar <- b
#         k <- 1L
#         # Update the VCV matrix
#         while (k <= (m + 1)) {
#             i <- seg[k] + 1
#             j <- seg[k + 1]
#             tmpy <- maty[((i - 1)*n+1):(j * n), , drop = FALSE]
#             tmpx <- matz[((i - 1)*n+1):(j * n), , drop = FALSE]
#             res <- tmpy - tmpx %*% bstar
#             umat <- reshapex(res, nrow(res)/n, n)
#             vvar <- crossprod(x = umat, y = NULL) / (j - i + 1)
#             vv[((k - 1)*n+1):(k * n), ] <- vvar
#             k <- k + 1L
#         }
        
#         k <- 1L
#         # Create the block diagonal matrix of the covariances
#         while (k <= (m + 1)) {
#             i <- seg[k] + 1
#             j <- seg[(k + 1)]
#             ibigv[((i - 1)*n+1):(j * n), ((i - 1)*n+1):(j * n)] <- eye(j - i + 1) %x% qr.solve(vv[((k - 1)*n+1):(k * n), ])
#             k <- k + 1L
#         }
        
#         b <- qr.solve(t(matz) %*% ibigv %*% matz, t(matz) %*% ibigv %*% maty)
#         itr <- itr + 1L
#         if (itr == max_iter) {
#             warning("Iteration has reached the upper bound.\n", call. = FALSE)
#             break
#         }
#     }
        
#     nbeta <- R %*% b
#     nvv <- vv
#     k <- 2L
#     while (k <= (m + 1)) {
#         nbeta <- rbind(nbeta, b)
#         k <- k + 1L
#     }
# }
    
# Changes in coefficients only
if ((brv == 1) && (brbeta != 1)) {
    k <- 1L
    vvar <- 0
    while (k <= (m + 1)) {
        i <- seg[k] + 1
        j <- seg[k + 1]
        tmpy <- maty[((i - 1)*n+1):(j * n), , drop = FALSE]
        tmpx <- matz[((i - 1)*n+1):(j * n), , drop = FALSE]
        res <- tmpy - tmpx %*% t(beta[k, , drop = FALSE])
        umat <- reshapex(res, nrow(res)/n, n)
        vvar <- vvar + crossprod(x = umat, y = NULL)
        k <- k + 1L
    }
    
    vvar <- vvar / bigt
    ibigv <- eye(bigt) %x% qr.solve(vvar)
    pmatz <- pzbar(matz, m, br, bigt)
    b <- qr.solve(t(pmatz) %*% ibigv %*% pmatz, t(pmatz) %*% ibigv %*% maty)
    # Now iterate until convergence is reached
    bstar <- b + 10
    itr <- 1L
    while ((maxc(abs(bstar - b)) > 1E-6) && (itr < max_iter)) {
        bstar <- b
        k <- 1L
        res <- maty - pmatz %*% bstar
        umat <- reshapex(res, nrow(res)/n, n)
        vvar <- crossprod(x = umat, y = NULL) / bigt
        ibigv <- eye(bigt) %x% qr.solve(vvar)
        b <- qr.solve(t(pmatz) %*% ibigv %*% pmatz, t(pmatz) %*% ibigv %*% maty)
        itr <- itr + 1
        if (itr == max_iter) {
            warning("Iteration has reached the upper bound.\n", call. = FALSE)
            break
        }
    }
    
    nbeta <- R %*% b
    nvv <- vvar
    k <- 2L
    while (k <= (m + 1)) {
        nvv <- rbind(nvv, vvar)
        k <- k + 1L
    }
}

# return estimated coefficients and covariances
return (list(nbeta = nbeta, nvv = nvv))

}

# Construct the diagonal partition of matz with m breaks at dates bb
pzbar <- function(matz, m, bb, bigt) {

    # Get some useful values
    nt <- nrow(matz)
    q1 <- ncol(matz)
    n <- nt / bigt
    
    if (m == 0) {
        stop("m == 0, no break is allowed.\n")
    } else {
        matzb <- matrix(data = 0, nrow = nt, ncol = (m + 1) * q1)
        matzb[1:(bb[1] * n), 1:q1] <- matz[1:(bb[1] * n), ]
        i <- 2L
        while (i <= m) {
            matzb[(bb[i - 1]*n+1):(bb[i] * n), ((i - 1)*q1+1):(i * q1)] <- matz[(bb[i - 1]*n+1):(bb[i] * n), ]
            i <- i + 1L
        }
        
        matzb[(bb[m]*n+1):nt, (m*q1+1):((m + 1) * q1)] <- matz[(bb[m]*n+1):nt, ]
    }
    return (matzb)
}    

# Estimate robust standard errors
correct <- function(res, prewhit) {
    
    # Get some useful values
    d <- ncol(res)
    nt <- nrow(res)
    bhat <- matrix(data = 0.0, nrow = d, ncol = 1)
    bmat <- matrix(data = 0.0, nrow = d, ncol = d)
    vstar <- matrix(data = 0.0, nrow = (nt - 1), ncol = d)
    vmat <- res
    
    # Apply pre-whitening?
    if (prewhit) {
        for (i in 1L:d) {
            bhat <- qr.coef(qr = qr(vmat[1:(nt - 1), , drop = FALSE]), y = vmat[2:nt, i, drop = FALSE])
            bmat[, i] <- bhat
            vstar[, i] <- vmat[2:nt, i, drop = FALSE] - vmat[1:(nt - 1), , drop = FALSE] %*% bhat
        }
    
        # Call the kernel on the residuals
        jh <- jhatpr(vstar)
        
        # Re-colour
        I_bhat_i <- qr.solve(eye(d) - bmat)
        hac <- I_bhat_i %*% jh %*% t(I_bhat_i)
    } else {
        hac <- jhatpr(vmat)
    }
    return (hac)
}

# Compute long-run VCV matrix of vmat
jhatpr <- function(vmat) {

    # Get a few useful values
    nt <- nrow(vmat)
    d <- ncol(vmat)
    jhat <- matrix(data = 0.0, nrow = d, ncol = d)
    
    # Call the automatic bandwidth selector
    st <- bandw(vmat)
    
    # HAC correction: Lag 0
    jhat <- crossprod(x = vmat, y = NULL) # vmat'vmat
    
    # HAC correction: Lag 1,...,T-1
    # Forward summation...
    for (j in 1L:(nt - 1)) { 
        jhat <- jhat + kern(j/st) * crossprod(x = vmat[(j + 1):nt, ], y = vmat[1:(nt - j), ])
    }
    # ...and Backward summation.
    for (j in 1L:(nt - 1)) {
        jhat <- jhat + kern(j/st) * crossprod(x = vmat[1:(nt - j), ], y = vmat[(j + 1):nt, ])
    }
    
    # Do small sample correction as specified by Andrews (1991)
    jhat <- jhat / (nt - d)
    return (jhat)
}

# Compute Andrews (1991) automatic bandwidth using and AR(1) estimated via OLS
bandw <- function(vhat) {

    # Get some useful values
    nt <- nrow(vhat)
    d <- ncol(vhat)
    a2n <- 0.0
    a2d <- 0.0
    
    for (j in 1L:d) {
        v <- vhat[, j, drop = FALSE]
        vt <- v[2:nt]
        vt_1 <- v[1:(nt - 1)]
        rho_j <- sum(vt_1 * vt) / sum(vt_1 * vt_1); # NB: We use ||x'y|| / ||x'x|| for speed
        et <- vt - (rho_j * vt_1);
        s2 <- sum(et * et) / (nt - 1);
        a2n <- a2n + ((4 * rho_j * rho_j * s2 * s2) / (1 - rho_j)^8)
        a2d <- a2d + ((s2 * s2) / (1 - rho_j)^4)
    }
    
    a2 <-a2n / a2d
    st <- 1.3221 * ((a2 * nt)^0.2)
    return (st)
}

# Compute Quadratic Spectral Kernel
kern <- function(x) {
    del <- (6.0 * pi * x) / 5.0
    ker <- 3.0 * (sin(del) / del - cos(del)) / (del * del)
    return (ker)
}

make_br_vec <- function(x)
{
    return (matrix(data = x, nrow = length(x)))
}

# Function to print out confidence sets information (length and included observation indices)
print_csm_info <- function(csm, level) {
    
    m <- ncol(csm)
    bigt <- nrow(csm)    
    
    cat("\n-------------------------------------------------------------------------\n")
    cat(sprintf("Confidence Set Level = %d%%\n\n", level))
    for(i in 1L:m) {
        csm_idx <- which(!is.na(csm[, i]), arr.ind = TRUE) 
        cat(sprintf("Break %d confidence set's length is %d observations and includes:\n", 
            i, length(csm_idx)))
        print.default(csm_idx)
        cat('\n')
    }

}

# Print method for class "invLR"
print.invLR <- function(x, digits = max(3L, getOption("digits") - 3L), print_all = TRUE, ...) {
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    
    cat("\n\tInverted LR Confidence Sets (Eo and Morley 2015):\n")
    
    cat("\nConfidence Sets:\n")
    print_csm_info(x$csm99, level = 99)
    
    print_csm_info(x$csm95, level = 95)
    
    print_csm_info(x$csm90, level = 90)
    
    if (print_all) {
        cat("\nCritical Values (kappa):\n")
        print.default(x$kmat)
        
        cat("\nRestricted Beta:\n")
        print.default(x$rbetast)
        
        cat("\nRestricted Variance:\n")
        print.default(x$rvvst)
    }
    cat('\n')
}

# Plot method for class "invLR"
plot.invLR <- function(x, col_y = "black", lwd_y = 1, lwd_cs = 3, lwd_br = 2, col_af = 0.1) {
    
    # Set up plotting device
    op <- par(mar = c(3, 3, 2, 1) + 0.1, oma = c(0, 0, 3, 0), bty = 'o', tcl = -0.5, 
        lty = 1, lwd = 1, xaxs = 'r', yaxs = 'r', mgp = c(3, 1, 0), 
        cex.main = 1.2, font.main = 2, cex.lab = 1.1, font.lab = 1, cex.axis = 1.1, font.axis = 1)
    
    on.exit(par(op))
    
    m <- ncol(x$csm99)
    y_names <- colnames(x$y)
    if(is.null(y_names)) {
        y_names <- paste0('Y', 1:ncol(x$y))
    }
    delta <- 1.5

    layout(matrix(1L:(3L * ncol(x$y)), nrow = 3L, ncol = ncol(x$y), byrow = FALSE))
    
    for (i in 1L:ncol(x$y)) {
        
        # Gets upper and lower bounds
        ub <- max(x$y[, i]) * delta
        lb <- min(x$y[, i]) * delta
    
        # 99% Confidence set
        plot(x$y[, i], type = 'n', col = NA, lty = 1, lwd = 1, 
            main = sprintf("99%% Level - %s", y_names[i]), ylab = "", xlab = "")
        
        # Confidence set markers and estimated break point
        for (j in 1L:m) {
            lines(ifelse(x$csm99[, j] == 1, ub, x$csm99[, j]), type = 'h', 
                col = adjustcolor(j+1, alpha.f = col_af), 
                lwd = lwd_cs, lend = "square", ylim = c(lb, ub))
            lines(ifelse(x$csm99[, j] == 1, lb, x$csm99[, j]), type = 'h',
                col = adjustcolor(j+1, alpha.f = col_af), 
                lwd = lwd_cs, lend = "square", ylim = c(lb, ub))
            abline(v = x$br[j], col = j+1, lty = 1, lwd = lwd_br)
        }
        
        abline(h = 0, col = "black", lty = 1, lwd = 1)
        lines(x$y[, i], type = 'l', col = col_y, lty = 1, lwd = lwd_y)
        box(which = "plot", col = "black", lty = 1, lwd = 1)
        
        # 95% Confidence set
        plot(x$y[, i], type = 'n', col = NA, lty = 1, lwd = 1, 
            main = sprintf("95%% Level - %s", y_names[i]), ylab = "", xlab = "")
        
        # Confidence set markers and estimated break point
        for (j in 1L:m) {
            lines(ifelse(x$csm95[, j] == 1, ub, x$csm95[, j]), type = 'h',
                col = adjustcolor(j+1, alpha.f = col_af),
                lwd = lwd_cs, lend = "square", ylim = c(lb, ub))
            lines(ifelse(x$csm95[, j] == 1, lb, x$csm95[, j]), type = 'h', 
                col = adjustcolor(j+1, alpha.f = col_af),
                lwd = lwd_cs, lend = "square", ylim = c(lb, ub))
            abline(v = x$br[j], col = j+1, lty = 1, lwd = lwd_br)
        }
        
        abline(h = 0, col = "black", lty = 1, lwd = 1)        
        lines(x$y[, i], type = 'l', col = col_y, lty = 1, lwd = lwd_y)
        box(which = "plot", col = "black", lty = 1, lwd = 1)
        
        # 90% Confidence set
        plot(x$y[, i], type = 'n', col = NA, lty = 1, lwd = 1, 
            main = sprintf("90%% Level - %s", y_names[i]), ylab = "", xlab = "")
        
        # Confidence set markers and estimated break point
        for (j in 1L:m) {
            lines(ifelse(x$csm90[, j] == 1, ub, x$csm90[, j]), type = 'h', 
                col = adjustcolor(j+1, alpha.f = col_af),
                lwd = lwd_cs, lend = "square", ylim = c(lb, ub))
            lines(ifelse(x$csm90[, j] == 1, lb, x$csm90[, j]), type = 'h',
                col = adjustcolor(j+1, alpha.f = col_af),
                lwd = lwd_cs, lend = "square", ylim = c(lb, ub))
            abline(v = x$br[j], col = j+1, lty = 1, lwd = lwd_br)
        }
        
        abline(h = 0, col = "black", lty = 1, lwd = 1)        
        lines(x$y[, i], type = 'l', col = col_y, lty = 1, lwd = lwd_y)
        box(which = "plot", col = "black", lty = 1, lwd = 1)
        
    }
    
    # Main plot title
    mtext("Inverted LR Confidence Sets", side = 3, outer = TRUE, cex = 1.2, font = 2)
}

# EOF








####################################################################################################
# Inverted Likelihood Ratio R Script
####################################################################################################
#
# This file provides the main function to calculate confidence sets for break dates as 
# documented in:
#
# Yunjong Eo and James Morley (2015)
#   "Likelihood-Ratio-Based Confidence Sets for the Timing of Structural Breaks",
#   Quantitative Economics, Volume 6, Issue 2, pages 463–-497, DOI: 10.3982/QE186.
#
# These codes are a direct translation of the original GAUSS procs written by:
# Yunjong Eo and James Morley (2015)
#
# R emulation of GAUSS procs, printing and plotting codes by Luke Hartigan, 28-07-2016
####################################################################################################

invLR_cs <- function(y, z, S, m, bigt, R, ...) {
  UseMethod("invLR_cs")
}

invLR_cs.default <- function(y, z, S, m, bigt, br, hetq, vauto, brv, brbeta, R, prewhit, trm) {

    # Get some useful quantities
    y <- as.matrix(y)
    z <- as.matrix(z)
    n <- ncol(y)
    
    # Set up regression matrices
    tmp <- transf(y, z, S)
    maty <- tmp$maty
    matz <- tmp$matz
    
    # Estimate the restricted model
    tmp <- restim(maty, matz, bigt, n, m, br, R, brbeta, brv)
    rbetast <- tmp$nbeta
    rvvst <- tmp$nvv
    likst <- estim_lik(maty, matz, n, m, br, rbetast, rvvst)
    
    # Compute scale factors
    wmat <- ws(maty, matz, m, n, bigt, br, rbetast, rvvst, hetq, vauto, brv, brbeta, prewhit)
    
    kmat <- matrix(data = 0.0, nrow = m, ncol = 3L)
    
    j <- 1L
    while (j <= m) {
        w1 <- wmat[j, 1L]
        w2 <- wmat[j, 2L]
        
        kmat[j, 1L] <- cv(w1, w2, 0.99)
        kmat[j, 2L] <- cv(w1, w2, 0.95)
        kmat[j, 3L] <- cv(w1, w2, 0.90)
        
        j <- j + 1L
    }
    
    csm99 <- csm95 <- csm90 <- matrix(data = NA_real_, nrow = bigt, ncol = m)
    
    gp <- round(trm * bigt)
    
    brm <- rbind(0, br, nrow(maty)/n)
    
    k0 <- 1L
    
    while (k0 <= m) {
        i <- brm[k0] + gp
        j <- brm[k0 + 2L] - gp
        
        cn <- i
        while (cn <= j) {
            brc <- br
            brc[k0] <- cn
            tmp <- restim(maty, matz, bigt, n, m, brc, R, brbeta, brv)
            nbeta <- tmp$nbeta
            nvv <- tmp$nvv
            lik <- estim_lik(maty, matz, n, m, brc, nbeta, nvv)
            
            if ((likst - lik) <= kmat[k0, 1L]) {
                csm99[cn, k0] <- 1L
            }
            
            if ((likst - lik) <= kmat[k0, 2L]) {
                csm95[cn, k0] <- 1L
            }
            
            if ((likst - lik) <= kmat[k0, 3L]) {
                csm90[cn, k0] <- 1L
            }
            
            cn <- cn + 1L
        }
        
        k0 <- k0 + 1L
    }
    
    # Collect results as a list object
    invlr <- list(y = y, br = br, csm99 = csm99, csm95 = csm95, csm90 = csm90, 
        kmat = kmat, rbetast = rbetast, rvvst = rvvst)
    
    # Return results with class "invLR"
    invlr$call <- match.call()
    class(invlr) <- "invLR"
    return(invlr)
}

# Function to calculate scale factors for the limit distribution
# See Eo and Morley (2015) for notations (B1, B2, Q1, Q2, Pi1, Pi2, Omega1, Omega2, Psi1, Psi2, Gam1, Gam2)
ws <- function(maty, matx, m, n, bigt, br, betax, vv, hetq, vauto, brv, brbeta, prewhit) {
    wmat <- matrix(data = 0.0, nrow = m, ncol = 2)
    diagx <- pzbar(matx, m, br, bigt)
    br <- rbind(br, (nrow(maty) / n))
    res <- maty - (diagx %*% betax)
    res <- reshapex(res, nrow(res) / n, n)
    
    # Get standardised residuals
    j <- 1L
    while (j <= (m + 1)) {
        if (j == 1L) {
            res[1:br[j], ] <- res[1:br[j], ] %*% qr.solve(sqrm(vv[1:(j * n),]))
        } else {
            res[(br[j - 1] + 1):br[j], ] <- res[(br[j - 1] + 1):br[j], ] %*% qr.solve(sqrm(vv[((j - 1) * n + 1):(j * n), ]))
        }
    j <- j + 1L
    }
    
    j <- 1L
    while (j <= m) {
        vvar1 <- vv[((j - 1) * n + 1):(j * n), ]
        vvar2 <- vv[(j * n + 1):((j + 1) * n), ]
        B1 <- sqrm(vvar1) %*% qr.solve(vvar2) %*% (vvar2 - vvar1) %*% qr.solve(sqrm(vvar1))
        B2 <- sqrm(vvar2) %*% qr.solve(vvar1) %*% (vvar2 - vvar1) %*% qr.solve(sqrm(vvar2))
        Q1 <- Q2 <- Pi1 <- Pi2 <- 0.0
        if (j == 1L) {
            i <- 1L
        } else {
            i <- br[j - 1] + 1L
        }
        if (hetq == 0) {
            # Change in variance and/or regression coefficients, 
            # allow for serial correlation and change in the distribution of the regressors
            if ((brv == 1) && (vauto == 1)) { 
                # Construct estimate of Pi1
                tempmat <- matrix(data = 0.0, nrow = (br[j] - i + 1), ncol = ncol(matx))
                k <- 1L
                while (k <= (br[j] - i + 1)) {
                    tempmat[k, ] <- t(t(matx[(((i-1+k-1)*n+1):(i-1+k)*n), , drop = FALSE]) %*% qr.solve(vvar2) %*% sqrm(vvar1) %*% t(res[i-1+k, , drop = FALSE]))
                    k <- k + 1L
                }
                Pi1 <- correct(tempmat, prewhit)        # with scaling
                
                # Construct estimate of Omega1
                tempmat <- matrix(data = 0.0, nrow = (br[j] - i + 1), n^2)
                k <- 1L
                while (k <= (br[j] - i + 1)) {
                    tempmat[k, ] <- t(vec(crossprod(x = res[(i - 1 + k), , drop = FALSE], y = NULL) - eye(n)))
                    k <- k + 1L
                }
                Omega1 <- correct(tempmat, prewhit)     # with scaling
                
                while (i <= br[j]) {
                    Q1 <- Q1 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                
                if (j == 1) {
                    Q1 <- Q1 / br[j]
                } else {
                    Q1 <- Q1 / (br[j] - br[j - 1])
                }
                
                i <- br[j] + 1L
                
                # Construct estimate of Pi2
                tempmat <- matrix(data = 0.0, nrow = (br[j + 1] - i + 1), ncol = ncol(matx))
                k <- 1L
                while (k <= (br[j + 1] - i + 1)) {
                    tempmat[k, ] <- t(t(matx[(((i-1+k-1)*n+1):(i-1+k)*n), , drop = FALSE]) %*% qr.solve(vvar1) %*% sqrm(vvar2) %*% t(res[i-1+k, , drop = FALSE]))
                    k <- k + 1L
                }
                Pi2 <-  correct(rempmat, prewhit)       # with scaling
                
                # Construct estimate of Omega2
                tempmat <- matrix(data = 0.0, nrow = (br[j + 1] - i + 1), n^2)
                k <- 1L
                while (k <= (br[j + 1] - i + 1)) {
                    tempmat[k, ] <- t(vec(crossprod(x = res[(i - 1 + k), , drop = FALSE], y = NULL) - eye(n)))
                    k <- k + 1L
                }
                Omega2 <- correct(tempmat, prewhit)     # with scaling
                
                while (i <= br[j + 1]) {
                    Q2 <- Q2 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar1) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                Q2 <- Q2 / (br[j + 1] - br[j])
                
                diffb <- betax[(j * ncol(matx) + 1):((j + 1) * ncol(matx))] - betax[((j - 1) * ncol(matx) + 1):(j * ncol(matx))]
                Psi1 <- 0.5 * sumc(diag(B1 %*% B1)) + t(diffb) %*% Q1 %*% diffb
                Psi2 <- 0.5 * sumc(diag(B2 %*% B2)) + t(diffb) %*% Q2 %*% diffb
                Gam1 <- sqrt(0.25 * t(vec(B1)) %*% Omega1 %*% vec(B1) + t(diffb) %*% Pi1 %*% diffb)
                Gam2 <- sqrt(0.25 * t(vec(B2)) %*% Omega2 %*% vec(B2) + t(diffb) %*% Pi2 %*% diffb)
            } 
            # Change in the variance and/or in the regression coefficients
            if ((brv == 1) && (vauto == 0)) {
                # Construct estimate of Omega1
                Omega1 <- 0
                while (i <= br[j]) {
                    Omega1 <- Omega1 + tcrossprod(x = vec(crossprod(x = res[i, , drop = FALSE], y = NULL) - eye(n)), y = NULL)
                    Q1 <- Q1 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    Pi1 <- Pi1 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% vvar1 %*% qr.solve(vvar2) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                
                if (j == 1L) {
                    Omega1 <- Omega1 / br[j]
                    Q1 <- Q1 / br[j]
                    Pi1 <- Pi1 / br[j]
                } else {
                    Omega1 <- Omega1 / (br[j] - br[j - 1])
                    Q1 <- Q1 / (br[j] - br[j - 1])
                    Pi1 <- Pi1 / (br[j] - br[j - 1])
                }
                
                i <- br[j] + 1L
                
                # Construct estimate of Omega2
                Omega2 <- 0
                while (i <= br[j + 1]) {
                    Omega2 <- Omega2 + tcrossprod(x = vec(crossprod(x = res[i, , drop = FALSE], y = NULL) - eye(n)), y = NULL)
                    Q2 <- Q2 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar1) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    Pi2 <- Pi2 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar1) %*% vvar2 %*% qr.solve(vvar1) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1
                }
                
                Omega2 <- Omega2 / (br[j] - br[j - 1])
                Q2 <- Q2 / (br[j] - br[j - 1])
                Pi2 <- Pi2 / (br[j] - br[j - 1])
                
                diffb <- betax[(j * ncol(matx) + 1):((j + 1) * ncol(matx))] - betax[((j - 1) * ncol(matx) + 1):(j * ncol(matx))]
                Psi1 <- 0.5 * sumc(diag(B1 %*% B1)) + t(diffb) %*% Q1 %*% diffb
                Psi2 <- 0.5 * sumc(diag(B2 %*% B2)) + t(diffb) %*% Q2 %*% diffb
                Gam1 <- sqrt(0.25 * t(vec(B1)) %*% Omega1 %*% vec(B1) + t(diffb) %*% Pi1 %*% diffb)
                Gam2 <- sqrt(0.25 * t(vec(B2)) %*% Omega2 %*% vec(B2) + t(diffb) %*% Pi2 %*% diffb)
            }
            # Change in the coefficients only
            if ((brv == 0) && (brbeta == 1) && (vauto == 1)) {
                # Construct estimate of Pi1
                tempmat <- matrix(data = 0.0, nrow = bigt, ncol = ncol(matx))
                k <- 1L
                while (k <= bigt) {
                    tempmat[k, ] <- t(t(matx[((k - 1)*n+1):(k * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% sqrm(vvar1) %*% t(res[k, , drop = FALSE]))
                    k <- k + 1L
                }
                Pi1 <- correct(tempmat, prewhit)    # with scaling
                k <- 1L
                while (k <= bigt) {
                    Q1 <- Q1 + t(matx[((k - 1)*n+1):(k * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% matx[((k - 1)*n+1):(k * n), , drop = FALSE]
                    k <- k + 1L
                }
                Q1 <- Q1 / bigt
                #Pi2 <- Pi1
                #Q2 <- Q1
                diffb <- betax[(j * ncol(matx) + 1):((j + 1) * ncol(matx))] - betax[((j - 1) * ncol(matx) + 1):(j * ncol(matx))]
                Psi1 <- t(diffb) %*% Q1 %*% diffb
                Psi2 <- Psi1 # t(diffb) %*% Q2 %*% diffb
                Gam1 <- sqrt(t(diffb) %*% Pi1 %*% diffb)
                Gam2 <- Gam1 # Gam1 <- sqrt(t(diffb) %*% Pi2 %*% diffb)
            }
            # Change in the coefficients only
            if ((brv == 0) && (brbeta == 1) && (vauto == 0)) {
                k <- 1L
                while (k <= bigt) {
                    Q1 <- Q1 + t(matx[((k - 1)*n+1):(k * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% matx[((k - 1)*n+1):(k * n), , drop = FALSE]
                    Pi1 + t(matx[((k - 1)*n+1):(k * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% vvar1 %*% qr.solve(vvar2) %*% matx[((k - 1)*n+1):(k * n), , drop = FALSE]
                    k <- k + 1L
                }
                Q1 <- Q1 / bigt
                Pi1 <- Pi1 / bigt
                #Pi2 <- Pi1
                #Q2 <- Q1
                diffb <- betax[(j * ncol(matx) + 1):((j + 1) * ncol(matx))] - betax[((j - 1) * ncol(matx) + 1):(j * ncol(matx))]
                Psi1 <- t(diffb) %*% Q1 %*% diffb
                Psi2 <- Psi1 # t(diffb) %*% Q2 %*% diffb
                Gam1 <- sqrt(t(diffb) %*% Pi1 %*% diffb)
                Gam2 <- Gam1 # Gam1 <- sqrt(t(diffb) %*% Pi2 %*% diffb)
            }
        } else if (hetq == 1) {
            # Allow the distribution of the regressors to be different, we only need to separately treat vauto == [1 || 0],
            # for a partial break model, the corresponding coefficients are zero
            if (vauto == 1) {
                # Construct estimate of Pi1
                tempmat <- matrix(data = 0.0, nrow = (br[j]- i + 1), ncol = ncol(matx))
                k <- 1L
                while (k <= (br[j] - i + 1)) {
                    tempmat <- t(t(matx[((i-1+k-1)*n+1):((i-1+k)*n), , drop = FALSE]) %*% qr.solve(vvar2) %*% sqrm(vvar1) %*% t(res[(i-1+k), , drop = FALSE]))
                    k <- k + 1L
                }
                Pi1 <- correct(tempmat, prewhit)    # with scaling
                
                # Construct estimate of Omega1
                tempmat <- matrix(data = 0.0, nrow = (br[j]-i+1), ncol = n^2)
                k <- 1L
                while (k <= (br[j] - i + 1)) {
                    tempmat[k, ] <- t(vec(crossprod(x = res[(i-1+k), , drop = FALSE], y = NULL) - eye(n)))
                    k <- k + 1L
                }
                Omega1 <- correct(tempmat, prewhit) # with scaling
                while (i <= br[j]) {
                    Q1 <- Q1 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                
                if (j == 1L) {
                    Q1 <- Q1 / br[j]
                } else {
                    Q1 <- Q1 / (br[j] - br[j - 1])
                }
                
                i <- br[j] + 1
                
                # Construct estimate of Pi2
                tempmat <- matrix(data = 0.0, nrow = (br[j + 1]-i+1), ncol = ncol(matx))
                k <- 1L
                while (k <= (br[j + 1] - i + 1)) {
                    tempmat[k, ] <- t(t(matx[((i-1+k-1)*n+1):((i-1+k)*n), , drop = FALSE]) %*% qr.solve(vvar1) %*% sqrm(vvar2) %*% t(res[(i-1+k), , drop = FALSE]))
                    k <- k + 1L
                }
                Pi2 <- correct(tempmat, prewhit)    # with scaling
                
                # Construct estimate of Omega2
                tempmat <- matrix(data = 0.0, nrow = (br[j+1]-i+1), ncol = n^2)
                k <- 1L
                while (k <= (br[j + 1] - i + 1)) {
                    tempmat[k, ] <- t(crossprod(x = res[(i-1+k), , drop = FALSE], y = NULL) - eye(n))
                    k <- k + 1L
                }
                Omega2 <- correct(tempmat, prewhit) # with scaling
                while (i <= br[j + 1]) {
                    Q2 <- Q2 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar1) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                Q2 <- Q2 / (br[j + 1] - br[j])
                diffb <- betax[(j * ncol(matx) + 1):((j + 1) * ncol(matx))] - betxa[((j - 1) * ncol(matx) + 1):(j * ncol(matx))]
                Psi1 <- 0.5 * sumc(diag(B1 %*% B1)) + t(diffb) %*% Q1 %*% diffb
                Psi2 <- 0.5 * sumc(diag(B2 %*% B2)) + t(diffb) %*% Q2 %*% diffb
                Gam1 <- sqrt(0.25 * t(vec(B1)) %*% Omega1 %*% vec(B1) + t(diffb) %*% Pi1 %*% diffb)
                Gam2 <- sqrt(0.25 * t(vec(B2)) %*% Omega2 %*% vec(B2) + t(diffb) %*% Pi2 %*% diffb)
            }
            # Change in the variance and/or in the regression coefficients
            if (vauto == 0) {
                # Construct estimate of Omega1
                Omega1 <- 0
                while (i <= br[j]) {
                    Omega1 <- Omega1 + tcrossprod(x = vec(crossprod(x = res[i, , drop = FALSE], y = NULL) - eye(n)), y = NULL)
                    Q1 <- Q1 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    Pi1 <- Pi1 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar2) %*% vvar1 %*% qr.solve(vvar2) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                
                if (j == 1) {
                    Omega1 <- Omega1 / br[j]
                    Q1 <- Q1 / br[j]
                    Pi1 <- Pi1 / br[j]
                } else {
                    Omega1 <- Omega1 / (br[j] - br[j - 1])
                    Q1 <- Q1 / (br[j] - br[j - 1])
                    Pi1 <- Pi1 / (br[j] - br[j - 1])
                }
                
                i <- br[j] + 1
                
                # Construct estimate of Omega2
                Omega2 <- 0
                while (i <= br[j + 1]) {
                    Omega2 <- Omega2 + tcrossprod(x = vec(crossprod(x = res[i, , drop = FALSE], y = NULL) - eye(n)), y = NULL)
                    Q2 <- Q2 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar1) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    Pi2 <- Pi2 + t(matx[((i - 1)*n+1):(i * n), , drop = FALSE]) %*% qr.solve(vvar1) %*% vvar2 %*% qr.solve(vvar1) %*% matx[((i - 1)*n+1):(i * n), , drop = FALSE]
                    i <- i + 1L
                }
                
                Omega2 <- Omega2 / (br[j + 1] - br[j])
                Q2 <- Q2 / (br[j + 1] - br[j])
                Pi2 <- Pi2 / (br[j + 1] - br[j])
                
                diffb <- betax[(j * ncol(matx) + 1):((j + 1) * ncol(matx))] - betax[((j - 1) * ncol(matx) + 1):(j * ncol(matx))]
                Psi1 <- 0.5 * sumc(diag(B1 %*% B1)) + t(diffb) %*% Q1 %*% diffb
                Psi2 <- 0.5 * sumc(diag(B2 %*% B2)) + t(diffb) %*% Q2 %*% diffb
                Gam1 <- sqrt(0.25 * t(vec(B1)) %*% Omega1 %*% vec(B1) + t(diffb) %*% Pi1 %*% diffb)
                Gam2 <- sqrt(0.25 * t(vec(B2)) %*% Omega2 %*% vec(B2) + t(diffb) %*% Pi2 %*% diffb)
            }
        }        
        
        w1 <- (Gam1^2) / Psi1
        w2 <- (Gam2^2) / Psi2
        wmat[j, ] <- cbind(w1, w2)   
        j <- j + 1L
    }
    
    # Return the scale factors
    return (wmat)
}

# Function to calculate the critical value for the inverted LR confidence set
cv <- function(w1, w2, lv, max_iter = 1000) {
    
    a <- 1 - lv
    w <- rbind(w1, w2)
    
    uw <- maxc(w)
    lw <- minc(w)
    
    uk <- -uw * log(1 - sqrt(lv))
    lk <- -lw * log(1 - sqrt(lv))
    
    k <- 0.5 * (uk + lk)  
    cr <- (1.0 - exp(-k / w1)) * (1.0 - exp(-k / w2))
    
    itr <- 1L
    while ((abs(lv - cr) >= 1.0E-5) && (itr <= max_iter)) {
        if ((cr - lv) > 0.0 ) {
            uk <- k
        } else {
            lk <- k
        }
        
        uc <- (1.0 - exp(-uk / w1)) * (1.0 - exp(-uk / w2))
        lc <- (1.0 - exp(-lk / w1)) * (1.0 - exp(-lk / w2))
        
        k <- ((lv - lc) * uk + (uc - lv) * lk) / (uc - lc)
        cr <- (1.0 - exp(-k / w1)) * (1.0 - exp(-k / w2))
        
        itr <- itr + 1L
        if (itr == max_iter) {
            warning("Iteration has reached the upper bound for finding cv of invLR.\n", call. = FALSE)
            break
        }
    }
    return (k)
}   

# Function to compute the log of the likelihood value given the break dates
estim_lik <- function(maty, matz, n, m, br, betax, vv) {
    
    qn <- ncol(matz)
    lik <- 0.0
    br <- rbind(0, br, (nrow(maty) / n))
    k <- 1L
    while (k <= m + 1L) {
        i <- br[k] + 1L
        j <- br[k + 1L]
        segx <- matz[((i - 1L) * n + 1L):(j * n), , drop = FALSE]
        segy <- maty[((i - 1L) * n + 1L):(j * n), , drop = FALSE]
        b <- betax[(qn * (k - 1L) + 1L):(qn * k), 1L]
        vvar <- vv[(n * (k - 1L) + 1L):(n * k), 1L:n]
        res <- segy - segx %*% b
        umat <- reshapex(res, nrow(res) / n, n)
        
        if (is.null(dim(vvar))) {
            lik <- lik - ((j - i + 1L) * n * 0.5 * log(2.0 * pi)) - (j - i + 1L)  * 0.5 * log(vvar) - 0.5 * 
                sumc(sumc(umat %*% qr.solve(vvar) * umat))
        } else {
            lik <- lik - ((j - i + 1L) * n * 0.5 * log(2.0 * pi)) - (j - i + 1L)  * 0.5 * log(det(vvar)) - 0.5 * 
                sumc(sumc(umat %*% qr.solve(vvar) * umat))
        }
        
        k <- k + 1L
    }
    return (lik)
}   