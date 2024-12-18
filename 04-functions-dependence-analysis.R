## FUNCTIONS FOR THE METHODOLOGY IN SECTION 1.3.3.1 (ASSESING AND ESTIMATING EXTREMAL DEPENDENCE)


mergetz <- Rcpp::cppFunction(
    'Rcpp::List mergetzf(Rcpp::NumericVector t, Rcpp::NumericVector z, double prec) {
 int i = 0;
 int j = 0;
 int k = 0;
 int n = t.size();
 double tmp = 0.0;
 NumericVector nt(n);
 NumericVector nz(n);
 	while ( i < n - 1 )
	{
		if ( t[i+1] - t[i] < prec ) {
			j++;
		} else {
			if ( j > 0 ){
				for(int l = (i - j + 1); l < i + 1; l++){
					tmp += z[l];
				}
				nt[k] = t[i];
				nz[k] = tmp/((double) j);

				tmp = 0.0;
				j = 0;
			}else{
				nt[k] = t[i];
				nz[k] = z[i];
			}
			k++;
		}
		i++;
	}

	if ( j > 0 ){
		for(int l = (i - j + 1); l < i + 1; l++){
			tmp += z[l];
		}
		nt[k] = t[i];
		nz[k] = tmp/((double) j);
	} else {
		nt[k] = t[n-1];
		nz[k] = z[n-1];
	}
	k++;
	int nn = k;
  return Rcpp::List::create(Named("nt") = nt,
                            Named("nz") = nz,
                            Named("nn") = nn);
}')

tieAverage = function(x, n = nrow(x), prec = 1e-10) {
  x = x[order(x[, 1]), ]
  res <- mergetz(t = as.double(x[, 1]),
                 z = as.double(x[, 2]),
                 prec = as.double(prec))
  return(cbind(res[["nt"]], res[["nz"]])[1:res[["nn"]], ])
}

## Pseudo-observations for the A-plot that account for censoring.
# CAREFUL: It is only implemented when one variable is censored (X),
# censoring status (status==1: observed, status==0 censored).
# The estimator of Stute (1993) is used, \hat F_1 and \hat F_2 are its margins.

GetPseudos <- function(X, Y, status) {
  #V <- rank(Y)/length(Y)
  km <- survfit(Surv(X, status) ~ 1)
  ind <-
    as.numeric(apply(
      cbind(X),
      1,
      FUN = function(x) {
        which(km$time == x)
      }
    ))
  U <- 1 - km$surv[ind]
  j <- length(km$time) #number of distinct values of $X$
  weight.proxy <-
    (c(1 - km$surv[1], km$surv[1:(j - 1)] - km$surv[2:j])) * length(X) / (km$n.event[1:j] +
                                                                            0.5 * (km$n.event[1:j] == 0))
  weights <- status * weight.proxy[ind]
  N <- length(Y)
  # The pseudos of Y are calculated as G_n(Y_i), where G_n is the margin of the bivariate KM estimator in Example 1 in Gribkova and Lopez (2015)
  V <- numeric(N)
  for (i in 1:N) {
    V[i] <- (1 / N) * (sum(weights * (Y <= Y[i])))
  }
  # V <- rank(Y)/length(Y) would be the "uncensored" way, the difference for the loss/ALAE data is minor
  return(list(U = U, V = V, weights = weights))
}

# This function calculates the weighted ECDF at the point (u,v)

Emp <- function(u, data, weights) {
  ind1 <- as.numeric(data[, 1] <= u[1])
  ind2 <- as.numeric(data[, 2] <= u[2])
  (1 / nrow(data)) * sum(ind1 * ind2 * weights)
}

# Calculating the points (T_i,Z_i) that constitute the A-plot.
# Censoring is again allowed in one variable only. This variable is stored in the first column of x, and its censoring status is "status".

eA <- function(x,
               w1 = 0,
               w2 = 0,
               n = nrow(x),
               status = rep(1, n)) {
  n <- nrow(x)
  out <- GetPseudos(x[, 1], x[, 2], status)
  u <- cbind(out$U, out$V)
  t <- numeric(n)
  A <- numeric(n)
  for (i in 1:n) {
    den <- u[i, 1] * u[i, 2]
    if (den == 1) {
      den <- den - 1 / (2 * n)
    }
    t[i] <- log(u[i, 2]) / log(den)
    A[i] <- log(Emp(u[i, ], data = u, weights = out$weights)) / log(den)
  }
  I <- apply(
    u,
    1,
    FUN = function(x) {
      (x[1] > w1) && (x[2] > w2)
    }
  )
  I2 = (u[, 1] * u[, 2] < 1)
  res = tieAverage(cbind(t[I & I2], A[I & I2]))
  return(cbind(res[-c(1, nrow(res)), ], nrow(res) - 2))
}

fitAspline <-
  function(x,
           lambda,
           norder = 4,
           n = nrow(x),
           w1 = 0,
           w2 = 0,
           status = rep(1, n),
           l = floor(sqrt(n)),
           tau = 0.5,
           tt = seq(0, 1, l = 1001))
  {
    y = eA(
      x,
      w1 = w1,
      w2 = w2,
      n = n,
      status = status
    )
    n = y[1, 3]
    knn = sort(c(0, y[, 1], 1))
    knn = unique(knn)
    kn = knn[seq(1, length(knn), l = l)]
    X = fda::bsplineS(y[, 1], kn, norder = norder)
    lkn = ncol(X)
    R = fda::bsplineS(c(0, 1), kn, norder = norder)
    R = rbind(R, -R)
    r = c(1, 1, -1, -1)
    tR = fda::bsplineS(c(0, 1), kn, norder = norder, nderiv = 1)
    R = rbind(R, tR[1, ], -tR[2, ], -tR[1, ], tR[2, ])
    r = c(r, c(-1, -1, 0, 0))
    R = rbind(R, fda::bsplineS(kn, kn, norder = norder, nderiv = 2))
    r = c(r, rep(0, l))
    res0 = quantreg::rq(
      y[, 2] ~ X - 1,
      method = "fnc",
      R = R,
      r = r,
      tau = tau
    )

    tR = fda::bsplineS(kn, kn, norder = norder, nderiv = norder - 2)
    tst = tR
    tR = tR[-1, ] - tR[-nrow(tR), ]
    #  tR = tR[-c(1,nrow(tR)),]
    tR = rbind(tR, -tR)
    R = rbind(R, tR)

    #  lkn = ncol(X)
    #  Y = y[,2]
    #  tY = Y-X[,1]-X[,lkn]
    #  tX = X[,-c(1,lkn)]
    #  tR = R[,-c(1,lkn)]

    #  tr = c(r,rep(-lambda,nrow(tR)))
    #  tr = c(r, -lambda * rep(diff(kn)[-c(1,length(kn)-2)],2))
    tr = c(r, -lambda * rep(diff(kn), 2))
    res1 = quantreg::rq(
      y[, 2] ~ X - 1,
      method = "fnc",
      R = R,
      r = tr,
      tau = tau
    )
    bS = fda::bsplineS(tt, kn, norder = norder)
    bS1 = fda::bsplineS(tt, kn, norder = norder, nderiv = 1)
    bS2 = fda::bsplineS(tt, kn, norder = norder, nderiv = 2)
    return(list(
      "unconstrained" = list(
        "model" = res0,
        "A" = bS %*% res0$coef,
        "dA" = bS1 %*% res0$coef,
        "d2A" = bS2 %*% res0$coef,
        kn
      ),
      "constrained" = list(
        "model" = res1,
        "A" = bS %*% res1$coef,
        "dA" = bS1 %*% res1$coef,
        "d2A" = bS2 %*% res1$coef,
        kn
      ),
      "knots" = kn
    ))
  }

fitAsplineCV <-
  function(x,
           lambdas,
           norder = 4,
           n = nrow(x),
           w1 = 0,
           w2 = 0,
           status = rep(1, n),
           l = floor(sqrt(n)),
           tau = 0.5)
  {
    y = eA(
      x,
      w1 = w1,
      w2 = w2,
      n = n,
      status = status
    )
    n = y[1, 3]
    cv = matrix(0, nrow = n, ncol = length(lambdas))
    for (i in 1:n) {
      knn = sort(c(0, y[-i, 1], 1))
      knn = unique(knn)
      kn = knn[seq(1, length(knn), l = l)]
      X = fda::bsplineS(y[-i, 1], kn, norder = norder)
      R = fda::bsplineS(c(0, 1), kn, norder = norder)
      R = rbind(R, -R)
      r = c(1, 1, -1, -1)
      tR = fda::bsplineS(c(0, 1), kn, norder = norder, nderiv = 1)
      R = rbind(R, tR[1, ], -tR[2, ], -tR[1, ], tR[2, ])
      r = c(r, c(-1, -1, 0, 0))
      R = rbind(R, fda::bsplineS(kn, kn, norder = norder, nderiv = 2))
      r = c(r, rep(0, l))
      tR = fda::bsplineS(kn, kn, norder = norder, nderiv = norder - 2)
      tR = tR[-1, ] - tR[-nrow(tR), ]
      tR = rbind(tR, -tR)
      R = rbind(R, tR)
      for (j in 1:length(lambdas)) {
        tr = c(r, -lambdas[j] * rep(diff(kn), 2))
        res = quantreg::rq(
          y[-i, 2] ~ X - 1,
          method = "fnc",
          R = R,
          r = tr,
          tau = tau
        )
        cv[i,j] = cv[i,j] + abs(fda::bsplineS(y[i, 1], kn, norder = norder) %*% res$coef - y[i, 2])
      }
    }
    # Compute the average cross validation and standard error
    mean_cv <- colSums(cv)
    se_cv <- apply(cv, 2, sd) / sqrt(n)
    attr(mean_cv, "std. err") <- se_cv
    return(mean_cv)
  }

fitAcobs <-
  function(x,
           tt,
           n = nrow(x),
           l = floor(sqrt(n)),
           k = floor(n / 2),
           w1 = 0,
           w2 = 0,
           status = rep(1, n))
  {
    y = eA(
      x,
      w1 = w1,
      w2 = w2,
      n = n,
      status = status
    )
    n = y[1, 3]

    kn = c(0, quantile(y[, 1], seq(0, 1, l = l))[-c(1, l)], 1)
    # Constraints
    C1 = c(0, 0, 1)
    C2 = c(0, 1, 1)
    g = seq(1 / k, (k - 1) / k, 1 / k)
    L = cbind(rep(1, k - 1), g, pmax(g, 1 - g))
    M = rbind(C1, C2, L)
    # Fitting A with quadratic median B-spline smoothing
    Ac = cobs::cobs(
      c(0, y[, 1], 1),
      c(1, y[, 2], 1),
      constraint = "convex",
      lambda = -1,
      nknots = l,
      pointwise = M,
      print.warn = FALSE,
      print.mesg = FALSE
    )
    return(predict(Ac, z = tt, interval = "none")[, 2])
  }

# Generates a sample from the EV copula with the fitted B-spline estimator of A

EVCSample <-
  function(x,
           lambda,
           norder = 3,
           w1 = 0,
           w2 = 0,
           n = nrow(x),
           status = rep(1, n),
           l = floor(sqrt(n)),
           tau = 0.5) {
    y = eA(
      x,
      w1 = w1,
      w2 = w2,
      n = n,
      status = status
    )
    n = y[1, 3]
    knn = sort(c(0, y[, 1], 1))
    knn = unique(knn)
    kn = knn[seq(1, length(knn), l = l)]
    X = fda::bsplineS(y[, 1], kn, norder = norder)
    lkn = ncol(X)
    R = fda::bsplineS(c(0, 1), kn, norder = norder)
    R = rbind(R, -R)
    r = c(1, 1, -1, -1)
    tR = fda::bsplineS(c(0, 1), kn, norder = norder, nderiv = 1)
    R = rbind(R, tR[1, ], -tR[2, ], -tR[1, ], tR[2, ])
    r = c(r, c(-1, -1, 0, 0))
    R = rbind(R, fda::bsplineS(kn, kn, norder = norder, nderiv = 2))
    r = c(r, rep(0, l))
    res0 = quantreg::rq(
      y[, 2] ~ X - 1,
      method = "fnc",
      R = R,
      r = r,
      tau = tau
    )
    tR = fda::bsplineS(kn, kn, norder = norder, nderiv = norder - 2)
    tst = tR
    tR = tR[-1, ] - tR[-nrow(tR), ]
    tR = rbind(tR, -tR)
    R = rbind(R, tR)
    tr = c(r, -lambda * rep(diff(kn), 2))
    res1 = quantreg::rq(
      y[, 2] ~ X - 1,
      method = "fnc",
      R = R,
      r = tr,
      tau = tau
    )
    Ader <- function(t) {
      bS1 = fda::bsplineS(t, kn, norder = norder, nderiv = 1)
      bS1 %*% res1$coef
    }
    Asec <- function(t) {
      bS2 = fda::bsplineS(t, kn, norder = norder, nderiv = 2)
      bS2 %*% res1$coef
    }
    Ahat <- function(t) {
      bS = fda::bsplineS(t, kn, norder = norder)
      bS %*% res1$coef
    }
    gZ <- function(z) {
      1 + (1 - 2 * z) * (Ader(z) / Ahat(z)) + z * (1 - z) * (Asec(z) * Ahat(z) -
                                                               Ader(z) ^ 2) / Ahat(z) ^ 2
    }
    w <- seq(from = 0, to = 1, by = 0.001)
    C <- max(gZ(w))
    ##generate uniform random variables. Accept/Reject
    Z <- NULL
    while (length(Z) < n) {
      U <- runif(1, 0, 1)
      V <- runif(1, 0, 1)
      if (U < gZ(V) / C) {
        Z <- c(Z, V)
      }
    }
    pZ <- function(z) {
      (z * (1 - z) * Asec(z)) / (Ahat(z) * gZ(z))
    }
    Ber <- rbinom(n, 1, pZ(Z))
    U2 <- runif(n, 0, 1)
    V2 <- runif(n, 0, 1)
    W <- U2 * Ber + U2 * V2 * (1 - Ber)
    Y1 <- W ^ {
      Z / Ahat(Z)
    }
    Y2 <- W ^ {
      (1 - Z) / Ahat(Z)
    }
    Sample <- cbind(Y1, Y2)
    return(Sample)
  }#end of function

## FUNCTIONS FOR THE METHODOLOGY IN SECTION 1.3.3.2 (FITTING A PARAMETRIC COPULA MODEL)

# Calculating pseudo-observations for the estimation procedure (Kaplan--Meier estimator is used for X, ECDF for Y)
# CAREFUL: It is only implemented when one variable is censored (X), status is its censoring status (status==1: observed, status==0 censored).

GetPseudosEst <- function(X, Y, status) {
  #V <- rank(Y)/length(Y)
  km <- survfit(Surv(X, status) ~ 1)
  ind <-
    as.numeric(apply(
      cbind(X),
      1,
      FUN = function(x) {
        which(km$time == x)
      }
    ))
  #U <- 1-km$surv[ind]
  tie <- km$n.event[ind] > 1
  j <- length(km$time) #number of distinct values of $X$
  weight.proxy <-
    (c(1 - km$surv[1], km$surv[1:(j - 1)] - km$surv[2:j])) * length(X) / (km$n.event[1:j] +
                                                                            0.5 * (km$n.event[1:j] == 0))
  weights <- status * weight.proxy[ind]
  N <- length(Y)
  # The pseudos of X are calculated as F_n(X_i), where F_n is the KM estimator, and the pseudos of Y use just the ecdf
  U <- numeric(N)
  U.left <- numeric(N)
  V <- numeric(N)
  V.left <- numeric(N)
  for (i in 1:N) {
    U[i] <- (1 / N) * (sum(weights * (X <= X[i])))
    U.left[i] <- (1 / N) * (sum(weights * (X < X[i])))
    V[i] <- (1 / N) * (sum((Y <= Y[i])))
    V.left[i] <- (1 / N) * (sum((Y < Y[i])))
  }
  # V <- rank(Y)/length(Y) would be the "uncensored" way, the difference for the loss/ALAE data is minor
  return(list(
    U = U,
    U.left = U.left,
    V = V,
    V.left = V.left,
    tie = tie
  ))
}

# Ery tie correction
#Implementing this only for the Gumbel Copula for simplicity

negloglik.cens.Ery <- function(theta, U, U.left, V, V.left, status) {
  cop <- gumbelCopula(theta) #alter this line for other families
  uncens <- cbind(U.left, U, V.left, V)[status == 1, ]
  cont.uncens <-
    pCopula(uncens[, c(2, 4)], cop) - pCopula(uncens[, c(1, 4)], cop) - pCopula(uncens[, c(2, 3)], cop) +
    pCopula(uncens[, c(1, 3)], cop)
  # careful: the diffs can be 0 when the observations are censored
  cens <- cbind(V, U)[status == 0, ]
  cont.cens <- 1 - cCopula(cens, cop)[, 2]
  - sum(c(log(cont.uncens), log(cont.cens)))
}

# Li et al. tie correction
#Implementing this only for the Gumbel Copula for simplicity

negloglik.cens.Li <- function(theta, U, U.left, V, V.left, status, tie) {
  cop <- gumbelCopula(theta) ##alter this line for other families
  n <- length(U)
  uncens.tie <- cbind(U.left, U, V*(n/(n+1)))[(status == 1) & (tie == 1), ]
  uncens.notie <- (cbind(U, V)[(status == 1) & (tie == 0), ]) * (n / (n +
                                                                        1))
  cont.uncens.notie <- dCopula(uncens.notie, cop)
  cont.uncens.tie <-
    cCopula(uncens.tie[, c(3, 2)], cop)[, 2] - cCopula(uncens.tie[, c(3, 1)], cop)[, 2]
  # careful: the diffs can be 0 when the observations are censored
  cens <- cbind(V, U)[status == 0, ]
  cont.cens <- 1 - cCopula(cens, cop)[, 2]
  - sum(c(
    log(cont.uncens.tie),
    log(cont.uncens.notie),
    log(cont.cens)
  ))
}

# Goodness of fit testing procedures

# Calculating the bootstrap sample in Section 1.3.3.1


sample.boot.n <-
  function(data,
           lambda,
           norder = 3,
           w1 = 0,
           w2 = 0,
           n = nrow(data),
           status = rep(1, n),
           l = floor(sqrt(n)),
           tau = 0.5,
           seed = 1,
           x,
           Fx,
           y,
           Fy,
           p) {
    set.seed(seed)
    U <-
      EVCSample(
        x = data,
        lambda = lambda,
        norder = norder,
        w1 = w1,
        w2 = w2,
        status = status,
        l = l,
        tau = tau
      )
    X.uncens.boot <-
      pmin(sapply(U[, 1], function(uval)
        min(x[Fx >= uval])), max(x) + 1) # the pmin is an adjustment for the case when the KM is not 1 at the largest observation (happens when the latter is censored)
    Y.boot <-
      pmin(sapply(U[, 2], function(uval)
        min(y[Fy >= uval])), max(y) + 1)#G$qfun(U[,2])
    perm <- unlist(sapply(Y.boot, function(yval)
      which(y == yval)))
    X.boot <- pmin(X.uncens.boot, p[perm])
    status.boot <- 1 - (X.boot < X.uncens.boot)
    return(cbind(X.boot, Y.boot, status.boot))
  }

# Calculating the test statistic in Section 1.3.3.1

Rn.n <-
  function(data,
           lambda,
           norder = 3,
           w1 = 0,
           w2 = 0,
           n = nrow(data),
           status = rep(1, n),
           l = floor(sqrt(n)),
           tau = 0.5) {
    tmp <- eA(
      x = data,
      status = status,
      w1 = w1,
      w2 = w2,
      n = n
    )
    m <-
      fitAspline(
        x = data,
        lambda = lambda,
        tt = tmp[, 1],
        norder = norder,
        n = n,
        w1 = w1,
        w2 = w2,
        l = l,
        status = status,
        tau = tau
      )
    sum((tmp[, 2] - m$con$A) ^ 2)
  }

# Calculating the bootstrap sample in Section 1.3.3.2

sample.boot <- function(n, theta, seed, x, Fx, y, Fy, p) {
  cop <- gumbelCopula(theta)
  set.seed(seed)
  U <- rCopula(n, cop)
  X.uncens.boot <-
    pmin(sapply(U[, 1], function(uval)
      min(x[Fx >= uval])), max(x) + 1) # the pmin is an adjustment for the case when the KM is not 1 at the largest observation (happens when the latter is censored)
  Y.boot <-
    pmin(sapply(U[, 2], function(uval)
      min(y[Fy >= uval])), max(y) + 1)#G$qfun(U[,2])
  perm <- unlist(sapply(Y.boot, function(yval)
    which(y == yval)))
  X.boot <- pmin(X.uncens.boot, p[perm])
  status.boot <- 1 - (X.boot < X.uncens.boot)
  return(cbind(X.boot, Y.boot, status.boot))
}

# Calculating the test statistic in Section 1.3.3.2

Rn <- function(x,
               w1 = 0,
               w2 = 0,
               n = nrow(x),
               status = rep(1, n)) {
  X <- x[, 1]
  Y <- x[, 2]
  pseudos <- GetPseudosEst(X, Y, status)
  theta.init <-
    summary(fitCopula(gumbelCopula(), pobs(cbind(
      pseudos$U, pseudos$V
    )), method = "mpl"))$coefficients[1]
  theta <-
    suppressWarnings(
      optim(
        theta.init,
        fn = negloglik.cens.Ery,
        U = pseudos$U,
        U.left = pseudos$U.left,
        V = pseudos$V,
        V.left = pseudos$V.left,
        status = status,
        method = "Brent",
        lower = 1,
        upper = 10
      )$par
    )
  tmp <- eA(
    x = x,
    status = status,
    w1 = w1,
    w2 = w2,
    n = n
  )
  sum((tmp[, 2] - A(gumbelCopula(theta), tmp[, 1])) ^ 2)
}

abvnonpar <- function (x = 0.5,
                       data,
                       epmar = TRUE,
                       method = c("cfg", "pickands"),
                       k = nrow(data) / 4,
                       ties.method = "random",
                       convex = FALSE,
                       rev = FALSE,
                       madj = 0,
                       kmar = NULL,
                       plot = FALSE,
                       add = FALSE,
                       lty = 1,
                       lwd = 1,
                       col = 1,
                       blty = 3,
                       blwd = 1,
                       xlim = c(0, 1),
                       ylim = c(0.5, 1),
                       xlab = "t",
                       ylab = "A(t)",
                       ...)
{
  if (mode(x) != "numeric" || any(x < 0, na.rm = TRUE) || any(x >
                                                              1, na.rm = TRUE))
    stop("invalid argument for `x'")
  method <- match.arg(method)
  epdata <-
    apply(data, 2, rank, na.last = "keep", ties.method = ties.method)
  nasm <- apply(data, 2, function(x)
    sum(!is.na(x)))
  epdata <- epdata / rep(nasm + 1, each = nrow(data))
  epdata <- -log(epdata)
  if (epmar)
    data <- epdata
  if (!epmar) {
    if (method == "pot") {
      if (any(k >= nasm))
        stop("k is too large")
      u1 <- sort(data[, 1], decreasing = TRUE)[k + 1]
      u2 <- sort(data[, 2], decreasing = TRUE)[k + 1]
      d1ab <- (data[, 1] > u1) & !is.na(data[, 1])
      d2ab <- (data[, 2] > u2) & !is.na(data[, 2])
      if (!is.null(kmar)) {
        data[d1ab, 1] <- mtransform(data[d1ab, 1], c(u1,
                                                     kmar))
        data[d2ab, 2] <- mtransform(data[d2ab, 2], c(u2,
                                                     kmar))
      }
      else {
        mle.m1 <- c(u1, fitted(fpot(data[d1ab, 1], threshold = u1)))
        mle.m2 <- c(u2, fitted(fpot(data[d2ab, 2], threshold = u2)))
        data[d1ab, 1] <- mtransform(data[d1ab, 1], mle.m1)
        data[d2ab, 2] <- mtransform(data[d2ab, 2], mle.m2)
      }
      data[d1ab, 1] <- -log(1 - k * data[d1ab, 1] / nasm[1])
      data[d2ab, 2] <- -log(1 - k * data[d2ab, 2] / nasm[2])
      data[!d1ab, 1] <- epdata[!d1ab, 1]
      data[!d2ab, 2] <- epdata[!d2ab, 2]
    }
    if (method != "pot") {
      if (!is.null(kmar)) {
        data <- mtransform(data, kmar)
      }
      else {
        if (!is.null(nsloc1)) {
          if (is.vector(nsloc1))
            nsloc1 <- data.frame(trend = nsloc1)
          if (nrow(nsloc1) != nrow(data))
            stop("`nsloc1' and data are not compatible")
          nslocmat1 <- cbind(1, as.matrix(nsloc1))
        }
        if (!is.null(nsloc2)) {
          if (is.vector(nsloc2))
            nsloc2 <- data.frame(trend = nsloc2)
          if (nrow(nsloc2) != nrow(data))
            stop("`nsloc2' and data are not compatible")
          nslocmat2 <- cbind(1, as.matrix(nsloc2))
        }
        mle.m1 <- fitted(fgev(data[, 1], nsloc = nsloc1,
                              std.err = FALSE))
        loc.mle.m1 <- mle.m1[grep("^loc", names(mle.m1))]
        if (is.null(nsloc1))
          loc.mle.m1 <- rep(loc.mle.m1, nrow(data))
        else
          loc.mle.m1 <- nslocmat1 %*% loc.mle.m1
        mle.m1 <- cbind(loc.mle.m1, mle.m1["scale"],
                        mle.m1["shape"])
        mle.m2 <- fitted(fgev(data[, 2], nsloc = nsloc2,
                              std.err = FALSE))
        loc.mle.m2 <- mle.m2[grep("^loc", names(mle.m2))]
        if (is.null(nsloc2))
          loc.mle.m2 <- rep(loc.mle.m2, nrow(data))
        else
          loc.mle.m2 <- nslocmat2 %*% loc.mle.m2
        mle.m2 <- cbind(loc.mle.m2, mle.m2["scale"],
                        mle.m2["shape"])
        data <- mtransform(data, list(mle.m1, mle.m2))
      }
    }
  }
  if (rev)
    data <- data[, 2:1]
  data <- na.omit(data)
  if (plot || add)
    x <- seq(0, 1, length = 100)
  d1 <- data[, 1]
  d2 <- data[, 2]
  sum1 <- sum(d1)
  slm1 <- sum(log(d1))
  sum2 <- sum(d2)
  slm2 <- sum(log(d2))
  nn <- nrow(data)
  nx <- length(x)
  mpmin <- function(a, b) {
    a[a > b] <- b[a > b]
    a
  }
  mpmax <- function(a, b) {
    a[a < b] <- b[a < b]
    a
  }
  if (method == "cfg") {
    if (!convex) {
      a <- numeric(nx)
      for (i in 1:nx)
        a[i] <- sum(log(mpmax((1 - x[i]) *
                                d1, x[i] * d2)))
      a <- (a - (1 - x) * slm1 - x * slm2) / nn
      a <- pmin(1, pmax(exp(a), x, 1 - x))
    }
    else {
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      for (i in 1:250)
        a[i] <- sum(log(mpmax((1 - x2[i]) *
                                d1, x2[i] * d2)))
      a <- (a - (1 - x2) * slm1 - x2 * slm2) / nn
      a <- pmin(1, pmax(exp(a), x2, 1 - x2))
      inch <- chull(x2, a)
      a <- a[inch]
      x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method = "linear")$y
    }
  }
  if (method == "pickands") {
    if (!convex) {
      a <- numeric(nx)
      if (madj == 2) {
        d1 <- d1 / mean(d1)
        d2 <- d2 / mean(d2)
      }
      for (i in 1:nx)
        a[i] <- sum(mpmin(d1 / x[i], d2 / (1 -
                                             x[i])))
      if (madj == 1)
        a <- a - x * sum1 - (1 - x) * sum2 + nn
      a <- nn / a
      a <- pmin(1, pmax(a, x, 1 - x))
    }
    else {
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      if (madj == 2) {
        d1 <- d1 / mean(d1)
        d2 <- d2 / mean(d2)
      }
      for (i in 1:250)
        a[i] <- sum(mpmin(d1 / x2[i], d2 / (1 -
                                              x2[i])))
      if (madj == 1)
        a <- a - x2 * sum1 - (1 - x2) * sum2 + nn
      a <- nn / a
      a <- pmin(1, pmax(a, x2, 1 - x2))
      inch <- chull(x2, a)
      a <- a[inch]
      x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method = "linear")$y
    }
  }
  if (method == "tdo") {
    if (!convex) {
      a <- numeric(nx)
      for (i in 1:nx)
        a[i] <- sum(mpmin(x[i] / (1 + nn *
                                    d1), (1 - x[i]) / (1 + nn * d2)))
      a <- 1 - a / (1 + log(nn))
      a <- pmin(1, pmax(a, x, 1 - x))
    }
    else {
      x2 <- seq(0, 1, length = 250)
      a <- numeric(250)
      for (i in 1:250)
        a[i] <- sum(mpmin(x2[i] / (1 + nn *
                                     d1), (1 - x2[i]) / (1 + nn * d2)))
      a <- 1 - a / (1 + log(nn))
      a <- pmin(1, pmax(a, x2, 1 - x2))
      inch <- chull(x2, a)
      a <- a[inch]
      x2 <- x2[inch]
      a <- approx(x2, a, xout = x, method = "linear")$y
    }
  }
  if (method == "pot") {
    a <- numeric(nx)
    rr <- rowSums(1 / data)
    rrk <- sort(rr, decreasing = TRUE)[k + 1]
    for (i in 1:nx)
      a[i] <- sum(mpmax(x[i] / (d1 * rr), (1 -
                                             x[i]) / (d2 * rr))[rr > rrk])
    a <- 2 / k * a
    a0 <- 2 / k * sum(1 / (d2 * rr)[rr > rrk])
    a1 <- 2 / k * sum(1 / (d1 * rr)[rr > rrk])
    a <- a + 1 - (1 - x) * a0 - x * a1
    a <- pmin(1, pmax(a, x, 1 - x))
  }
  if (plot || add) {
    if (!add) {
      plot(
        x,
        a,
        type = "n",
        xlab = xlab,
        ylab = ylab,
        xlim = xlim,
        ylim = ylim,
        ...
      )
      polygon(c(0, 0.5, 1), c(1, 0.5, 1), lty = blty, lwd = blwd)
    }
    lines(x,
          a,
          lty = lty,
          lwd = lwd,
          col = col)
    return(invisible(list(x = x, y = a)))
  }
  a
}
