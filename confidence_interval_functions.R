# functions to load

euc_norm <- function(x){sqrt(sum(x^2))}

fit_break_ols <- function(y, x, breaks) {
  N <- length(y)
  k <- ncol(x)
  breaks <- sort(unique(breaks))

  if (any(breaks <= 1 | breaks >= N)) {
    stop("Break points must be between 2 and N-1.")
  }

  seg_id <- cut(1:N, breaks = c(0, breaks, N), labels = FALSE, right = TRUE)

  X_big <- NULL
  for (s in 1:(length(breaks) + 1)) {
    Xs <- x
    Xs[seg_id != s, ] <- 0 
    X_big <- cbind(X_big, Xs)
  }

  coef <- solve(t(X_big) %*% X_big, t(X_big) %*% y)
  yhat <- X_big %*% coef
  resid <- y - yhat

  return(resid)
}


## utility functions for test statistics
## ols estimator
ols_e <- function(dv,indv){
    dv = as.matrix(dv)
    indv = as.matrix(indv)
    beta_hat = solve(t(indv)%*%indv) %*% (t(indv)%*%dv)
    return(beta_hat)
}

## residuals
resid <- function(dv,indv){
    N = length(dv)
    beta_hat = ols_e(dv,indv)
    residual = c()
    for(i in 1:N){
       residual[i] = dv[i]-as.matrix(indv)[i,]%*%beta_hat
    }
    return(residual)
}

confi_21 <- function(dv,indv,kappa){
    euc_norm <- function(x){sqrt(sum(x^2))}
    eps = resid(dv,indv)
    N = dim(indv)[1]
    d = dim(indv)[2]
    zseq = as.matrix(apply(indv, 2, function(x){x*eps}))
    zn = matrix(0, N, 1)
    for (j in 2:(N-1)){
       weight = 1/(j*(N-j))^kappa
       if (ncol(zseq)==1) {
         zn[j] = weight*euc_norm( sum(zseq[1:j,]) )
       }
       else {
         zn[j] = weight*euc_norm(colSums(zseq[1:j,]))
       }
    } 
    khat = which(zn==max(zn))
    return(khat)
}


matrixA <- function(indv){
  N = dim(indv)[1]
  d = dim(indv)[2]
  matA = t(as.matrix(indv))%*%as.matrix(indv) / N
  return(matA)
}




# check whether sqrt(var2) or var2, sqrt(var1) or var1
lim_sim_hete <- function(kappa, theta, N, nrep, a_minus, a_plus){
   nloc = 10000
   res1 = 1:nrep
   
   mkf1 = c(0, rep(((1-kappa)*(1-theta)+kappa*theta),nloc))
   mkf2 = c(0, rep(((1-kappa)*theta+kappa*(1-theta)),nloc))
   times = seq(0,N,length=(nloc+1))

   for(j in 1:nrep){
    x1 = (2*a_minus/(a_minus+a_plus))*BM(0,0,N,nloc) - abs(times)*mkf1
    
    x2 = (2*a_plus/(a_minus+a_plus) )*BM(0,0,N,nloc) - abs(times)*mkf2

    if(max(x1)>max(x2)){
      res1[j] = N*which(x1 == max(x1))/nloc
    } else{
      res1[j] = - N*which(x2 == max(x2))/nloc
    }
   }
   return( res1 )
}


# modified homo
sigma_istar <- function(dv,indv, khat, deltahat){
    euc_norm <- function(x){sqrt(sum(x^2))}
    N = dim(indv)[1]
    d = dim(indv)[2]
    
    eps = resid(dv,indv)
    # eps1 = resid(dv[1:khat],indv[1:khat,])
    # eps2 = resid(dv[(khat+1):N],indv[(khat+1):N,])
    # eps  = c(eps1,eps2)

    zseq = as.matrix(apply(indv, 2, function(x){x*eps}))
    Amat = matrixA(indv)

    dm = as.matrix(zseq)
    Dmat = (t(dm)%*%dm)/khat 
    
    ds = deltahat/euc_norm(deltahat)
    a1 = sqrt( t(ds)%*%Amat%*%Dmat%*%Amat%*%ds )
    a2 = a1

    sigmastar = ((a1+a2)/2 )^2 * (1/euc_norm(Amat%*%ds)^4)

    return(sigmastar)
}


# abrupt change standardization
sigma_step_correct <- function(dv,indv,khat,deltahat){
    
    eps = resid(dv,indv)
    # eps1 = resid(dv[1:khat],indv[1:khat,])
    # eps2 = resid(dv[(khat+1):length(dv)],indv[(khat+1):length(dv),])
    # eps  = c(eps1,eps2)

    N = dim(indv)[1]
    d = dim(indv)[2]
    zseq = as.matrix(apply(indv, 2, function(x){x*eps}))
    Amat = matrixA(indv)

    # dmd = as.matrix(apply(zseq,2,function(x){x-mean(x)}))
    Dmat = array(0, c(d, d, 2))
    dm = as.matrix(zseq[1:khat,])
    V = (t(dm)%*%dm)/khat 
    Dmat[,,1] = V
    dm = as.matrix(zseq[(khat+1):N,])
    V = (t(dm)%*%dm)/(N-khat)  
    Dmat[,,2] = V

    ds = deltahat/euc_norm(deltahat)

    sigma2 = array(0, 2)
    sigma2[1] = sqrt( t(ds)%*%Amat%*%Dmat[,,1]%*%Amat%*%ds )
    sigma2[2] = sqrt( t(ds)%*%Amat%*%Dmat[,,2]%*%Amat%*%ds )
    return(sigma2)
}

# smooth change standardization
sigma_step_smooth <- function(dv,indv,khat,deltahat, h_band){
    N = dim(indv)[1]
    d = dim(indv)[2]

    Amat = matrixA(indv)
    
    eps = resid(dv,indv)
    # eps1 = resid(dv[1:khat],indv[1:khat,])
    # eps2 = resid(dv[(khat+1):length(dv)],indv[(khat+1):length(dv),])
    # eps  = c(eps1,eps2)
    
    zseq = as.matrix(apply(indv, 2, function(x){x*eps}))
    
    Dmat = array(0, c(d, d, 2))
    left_start = khat-h_band
    if(left_start<0){left_start=1} # can be negative if khat<h_band
    dm = as.matrix(zseq[left_start:khat,])
    Dk = (t(dm)%*%dm)
    DV = Dk/ length(left_start:khat)
    Dmat[,,1] = DV
    
    right_end = khat+h_band
    if(right_end>N){right_end=N} # can be negative if khat<h_band
    dm = as.matrix(zseq[(khat+1):right_end,])
    Dk = (t(dm)%*%dm)
    DV = Dk/ length((khat+1):right_end)
    Dmat[,,2] = DV

    ds = deltahat/euc_norm(deltahat)

    sigma2 = array(0, 2)
    sigma2[1] = sqrt( t(ds)%*%Amat%*%Dmat[,,1]%*%Amat%*%ds ) # take square root
    sigma2[2] = sqrt( t(ds)%*%Amat%*%Dmat[,,2]%*%Amat%*%ds )
    return(sigma2)
}

sigma_standard <- function(dv,indv, khat, deltahat, g_minus_a, g_plus_a){    
    ds = deltahat/euc_norm(deltahat)
    Amat = matrixA(indv)

    sigmastar = ((g_minus_a+g_plus_a)/2 )^2 * (1/euc_norm(Amat%*%ds)^4)
    return(sigmastar)
}
 


ilr2015 <- function(y1, x1){
  bp.yx <- breakpoints(y1 ~ x1 )
 
  if (any(is.na(bp.yx$breakpoints)) == TRUE){
    inc = NA
    k_hat = NA
  }
  else{
    y = as.matrix(y1)
    if (ncol(x1)==1){
      z <- as.matrix( cbind(rep(1,length(x1)),x1) )
    }else{
      z <- as.matrix(x1)
    }

    # Specify some function options
    S <- eye(2)
    m <- length(bp.yx$breakpoints)
    bigt <- nrow(y)
    R <- eye(ncol(S) * (m + 1))
    trm <- 0.05
    brv <- 1
    brbeta <- 99
    vauto <- 0
    prewhit <- 0
    hetq <- 1

    br <- c( bp.yx[1]$breakpoints )
    br <- make_br_vec(br) # NB: 'br' must be a column vector for multiple breaks
    
    lh_ci = matrix(NA,2,m)
    # Run the code for the specified model and options provided
    results <- invLR_cs( y, z, S, m, bigt, br, hetq, vauto, brv, brbeta, R, prewhit, trm)
    for ( i in 1:m){
      lh_ci[1,i] = min(which(results$csm95[,i]==1))
      lh_ci[2,i] = max(which(results$csm95[,i]==1))
    }
 
    k_hat = br
  }

  return( list(lh_ci, k_hat ) )
}


 

# Fast residual sum of squares for y ~ X (no copies, QR via .lm.fit)
rss_ols <- function(y, X) {
  X <- as.matrix(X)
  fit <- .lm.fit(x = X, y = y)
  sum(fit$residuals^2)
}

# Return list of [start,end] (1-based, inclusive) given interior breaks (local)
segment_bounds_local <- function(n, breaks_local) {
  if (length(breaks_local)) {
    bk <- sort(unique(breaks_local))
    if (any(bk <= 1 | bk >= n)) stop("Local breaks must be in 2..(n-1).")
    starts <- c(1, bk + 1)
    ends   <- c(bk, n)
  } else {
    starts <- 1
    ends   <- n
  }
  data.frame(start = starts, end = ends)
}

# SSE for a (sub)sample y, X, split by *local* breaks
sse_fun <- function(y, X, breaks_local = integer(0)) {
  n <- length(y)
  if (n != nrow(X)) stop("y and X must have same number of rows.")
  segs <- segment_bounds_local(n, breaks_local)
  total <- 0
  for (i in seq_len(nrow(segs))) {
    a <- segs$start[i]; b <- segs$end[i]
    total <- total + rss_ols(y[a:b], X[a:b, , drop = FALSE])
  }
  total
}

## -------- optional: residuals for a *given* set of GLOBAL breaks --------
# (Only needed if you truly want residuals; SSE above is enough for selection.)
fit_break_residuals <- function(y, X, breaks_global = integer(0)) {
  N <- length(y)
  if (N != nrow(X)) stop("y and X must have same number of rows.")
  if (length(breaks_global)) {
    bk <- sort(unique(breaks_global))
    if (any(bk <= 1 | bk >= N)) stop("Global breaks must be in 2..(N-1).")
  } else {
    bk <- integer(0)
  }
  segs <- segment_bounds_local(N, bk)
  resid <- numeric(N)
  for (i in seq_len(nrow(segs))) {
    a <- segs$start[i]; b <- segs$end[i]
    fit <- .lm.fit(x = as.matrix(X[a:b, , drop = FALSE]), y = y[a:b])
    resid[a:b] <- fit$residuals
  }
  resid
}



greedy_breaks <- function(y, X, kappa,
                          max_breaks = 10,
                          min_seg_n = 20,
                          find_break = confi_21) {
  N <- length(y)
  if (N != nrow(X)) stop("y and X must have same number of rows.")
  segs <- data.frame(start = 1L, end = N)
  seg_sse <- rss_ols(y, X)

  breaks_global   <- integer(0)   # sorted set (bookkeeping)
  breaks_sequence <- integer(0)   # discovery order

  # optional: detailed log for each accepted split
  breaks_info <- data.frame(
    iter        = integer(0),
    k_glob      = integer(0),
    k_loc       = integer(0),
    start       = integer(0),
    end         = integer(0),
    gain        = numeric(0),
    sse_unsplit = numeric(0),
    sse_split   = numeric(0)
  )

  for (iter in seq_len(max_breaks)) {
    if (nrow(segs) == 0L) break

    proposals <- vector("list", nrow(segs))
    gains     <- rep(-Inf, nrow(segs))

    for (i in seq_len(nrow(segs))) {
      a <- segs$start[i]; b <- segs$end[i]
      n_i <- b - a + 1L
      if (n_i < 2L * min_seg_n + 1L) next

      y_i <- y[a:b]
      X_i <- X[a:b, , drop = FALSE]

      k_loc <- suppressWarnings(find_break(y_i, X_i, kappa))
      if (!is.numeric(k_loc) || length(k_loc) != 1L || is.na(k_loc)) next
      if (k_loc <= min_seg_n || (n_i - k_loc) <= min_seg_n) next

      sse_unsplit <- if (length(seg_sse) >= i && !is.na(seg_sse[i])) seg_sse[i] else rss_ols(y_i, X_i)
      sse_split   <- sse_fun(y_i, X_i, breaks_local = k_loc)  # LOCAL break
      gain <- sse_unsplit - sse_split

      proposals[[i]] <- list(k_loc = k_loc, a = a, b = b,
                             sse_split = sse_split, sse_unsplit = sse_unsplit)
      gains[i] <- gain
    }

    best_i <- which.max(gains)
    if (length(best_i) == 0L || !is.finite(gains[best_i]) || gains[best_i] <= 0) break

    pr <- proposals[[best_i]]
    a <- pr$a; b <- pr$b; k_loc <- pr$k_loc
    k_glob <- a + k_loc - 1L  # LOCAL -> GLOBAL

    # record in discovery order
    breaks_sequence <- c(breaks_sequence, k_glob)

    # insert into sorted set
    breaks_global <- sort(c(breaks_global, k_glob))

    # log details
    breaks_info <- rbind(breaks_info, data.frame(
      iter        = iter,
      k_glob      = k_glob,
      k_loc       = k_loc,
      start       = a,
      end         = b,
      gain        = gains[best_i],
      sse_unsplit = pr$sse_unsplit,
      sse_split   = pr$sse_split
    ))

    # replace chosen segment with its children
    left  <- data.frame(start = a,          end = k_glob)
    right <- data.frame(start = k_glob + 1, end = b)

    if (nrow(segs) == 1L) {
      segs <- rbind(left, right)
      seg_sse <- c(
        rss_ols(y[left$start:left$end],   X[left$start:left$end, , drop = FALSE]),
        rss_ols(y[right$start:right$end], X[right$start:right$end, , drop = FALSE])
      )
    } else {
      segs <- rbind(
        if (best_i > 1) segs[1:(best_i - 1), , drop = FALSE],
        left,
        right,
        if (best_i < nrow(segs)) segs[(best_i + 1):nrow(segs), , drop = FALSE]
      )
      seg_sse <- c(
        if (best_i > 1) seg_sse[1:(best_i - 1)] else NULL,
        rss_ols(y[left$start:left$end],   X[left$start:left$end, , drop = FALSE]),
        rss_ols(y[right$start:right$end], X[right$start:right$end, , drop = FALSE]),
        if (best_i < length(seg_sse)) seg_sse[(best_i + 1):length(seg_sse)] else NULL
      )
    }
  }

  # For convenience: map each sorted break to the iteration it was found
  breaks_sorted <- breaks_global
  breaks_order  <- match(breaks_sorted, breaks_sequence)  # 1 = found first, etc.

  list(
    breaks        = breaks_sorted,     # sorted (as before)
    breaks_sequence = breaks_sequence, # discovery order you asked for
    breaks_order  = breaks_order,      # per 'breaks', which iteration discovered it
    sse_total     = sse_fun(y, X, breaks_local = breaks_sorted),
    segments      = segs,
    breaks_info   = breaks_info        # detailed log (optional but useful)
  )
}
 

combinations_with_anchor <- function(breaks, anchor_idx = 1L, min_k = 1L, max_k = NULL) {
  n <- length(breaks)
  if (is.null(max_k)) max_k <- n
  if (anchor_idx < 1L || anchor_idx > n) stop("anchor_idx out of range.")
  if (min_k < 1L || max_k > n || min_k > max_k) stop("invalid min_k/max_k.")
  
  rest_idx <- setdiff(seq_len(n), anchor_idx)
  out <- list(); ctr <- 0L
  
  for (k in seq.int(min_k, max_k)) {     # total size including anchor
    if (k == 1L) {
      ctr <- ctr + 1L
      out[[ctr]] <- breaks[anchor_idx]
    } else {
      s <- k - 1L                        # how many from the rest
      if (length(rest_idx) >= s && s > 0L) {
        cmb <- combn(rest_idx, s, simplify = FALSE)
        for (idxs in cmb) {
          id <- sort(c(anchor_idx, idxs))      # keep original order
          ctr <- ctr + 1L
          out[[ctr]] <- breaks[id]
        }
      }
    }
  }
  out
}

info_criteria <- function(mstat, nchange, N){
  Kp = nchange *2 +1

  BIC_m = log(mstat) + 2*Kp * log(N)/N 
  # BIC_m = log(mstat) + Kp * log(N)/N
  return( BIC_m )
}


break_selection <- function(y,x,breaks){
  n_model = length(breaks)
  bic_m = matrix(NA,n_model,1)

  for (i in 1:n_model){
    resids = fit_break_ols(y, x, breaks[[i]])
    nchange = length(breaks[[i]])
    mstat = sum((resids)^2)
    N = length(resids)
    icr = info_criteria(mstat, nchange,  N)
    bic_m[i] = icr
  }

  bic_select = which(bic_m == min(bic_m))
  select_breaks = breaks[[bic_select]]
  return(sort(select_breaks))
}


# bre_out = greedy_breaks(y, x, kappa=1/2, max_breaks = 5, min_seg_n = 20)
# res_breaks <- combinations_with_anchor(bre_out$breaks_sequence)

# break_selection(y,x,res_breaks)


change_point_refinement <- function(y,x,breaks, kappa){
  split_by_breaks <- function(y, X, breaks) {
    if (length(y) != nrow(X)) stop("y and X must have the same number of rows.")
    N <- length(y)
    b <- sort(unique(as.integer(breaks)))
    if (length(b) < 2L) stop("Provide at least two breaks.")
    if (any(b <= 1L | b >= N)) stop("Breaks must be in 2..(N-1).")

    B <- c(1L, b, N)  # boundaries: 1, breaks..., N
    idx_pairs <- list()

    # sliding triples: (B[i] -> B[i+2]) for i = 1..(length(B)-3)
    if (length(B) >= 4L) {
      for (i in seq_len(length(B) - 3L)) {
        idx_pairs[[length(idx_pairs) + 1L]] <- c(B[i], B[i + 2L])
      }
    }
    #   tail windows to N
    idx_pairs[[length(idx_pairs) + 1L]] <- c(B[length(B) - 2L], B[length(B)])
    idx_pairs[[length(idx_pairs) + 1L]] <- c(B[length(B) - 1L], B[length(B)])

    # build outputs
    starts <- vapply(idx_pairs, `[[`, integer(1), 1L)
    ends   <- vapply(idx_pairs, `[[`, integer(1), 2L)
    ys <- mapply(function(a, b) y[a:b], starts, ends, SIMPLIFY = FALSE)
    Xs <- mapply(function(a, b) X[a:b, , drop = FALSE], starts, ends, SIMPLIFY = FALSE)

    names(ys) <- names(Xs) <- paste0("seg", seq_along(ys))
    list(bounds = data.frame(start = starts, end = ends), y = ys, x = Xs)
  }
  
  sample_split=split_by_breaks(y,x,breaks)
  
  r = length(breaks)
  refine_br = c()
  for (i in 1:(r+1)){
    refine_br[i] = confi_21(sample_split$y[[i]], sample_split$x[[i]], kappa)
  }
 
  refine_br = sample_split$bounds[,1]+refine_br-1
  
  if (length(refine_br) >= 2 && abs(refine_br[length(refine_br)] - refine_br[length(refine_br)-1]) < 20) {
    refine_br <- refine_br[-length(refine_br)]
  }

  return(refine_br)
}


 
ksc_confidence_interval <- function(y, x, refined_breaks, band_width, kappa = 1/2, alpha = 0.95){ 
  # alpha - nominal confidence coverage
  # kappa - weight parameter
  # band_width - long run variance estimator, e.g., ceiling( sqrt(sample_size) )
  split_by_breaks <- function(y, X, breaks) {
    N <- length(y)
    if (N != nrow(X)) stop("y and X must have same number of rows.")
    b <- sort(unique(as.integer(breaks)))
    m <- length(b)
    if (m < 2L) stop("Need at least two breaks to form these overlaps.")
    if (any(b <= 1L | b >= N)) stop("All breaks must be in 2..N-1.")

    starts <- c(1L, b[-m])      # 1, b1, b2, ..., b_{m-1}
    ends   <- c(b[-1L], N)      # b2, b3, ..., b_m, N

    idx_list <- Map(seq.int, starts, ends)
    y_list <- lapply(idx_list, function(ii) y[ii])
    X_list <- lapply(idx_list, function(ii) X[ii, , drop = FALSE])

    nm <- sprintf("[%d:%d]", starts, ends)
    names(y_list) <- names(X_list) <- nm

    list(starts = starts, ends = ends, idx = idx_list, y = y_list, X = X_list)
  }

  br_idx = split_by_breaks(y, x, refined_breaks)
  confi_interval = list()

  for(i in 1:length(refined_breaks)){ 
    y_sub = y[br_idx$starts[i]:br_idx$ends[i]]
    x_sub = as.matrix(x[br_idx$starts[i]:br_idx$ends[i],])
    N_sub = length(y_sub)
    khat_sub = refined_breaks[i] - br_idx$starts[i] + 1
    deltahat = ols_e(y_sub[1:khat_sub],x_sub[1:khat_sub,]) - ols_e(y_sub[(1+khat_sub):N_sub],x_sub[(1+khat_sub):N_sub,])
    varl=sigma_step_smooth(y_sub ,x_sub , khat_sub, deltahat, band_width)
    a_minus = varl[1]
    a_plus = varl[2]
    lim = lim_sim_hete(kappa=kappa, theta=khat_sub/N_sub, N=20000, nrep=5000, a_minus, a_plus )
    sigmastar = sigma_standard(y_sub,x_sub, khat_sub, deltahat, a_minus, a_plus )
    
    quan_lim = quantile(lim, c((1-alpha)/2, 1-(1-alpha)/2))
    rev_eng  = quan_lim*sigmastar/( euc_norm(deltahat)^2 )
    rev_eng[2] = max(0, rev_eng[2])
    confi = khat_sub - rev_eng
    confi_interval[[i]] = round(sort(confi) + br_idx$starts[i] - 1)
  }

  return(confi_interval)
}

