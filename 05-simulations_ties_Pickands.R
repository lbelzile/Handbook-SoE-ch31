setwd(this.path::here())
library(evd)
library(copula)
library(survival)
source("04-functions-dependence-analysis.R")
# Doing a small illustration to show that the uncorrected CFG and Pickands estimators
# are not unbiased when there are ties, but the B-spline estimator works well.

set.seed(1728)
n <- 5000
W <- rCopula(n, gumbelCopula(3))
Z <- cbind(W[, 1], qpois(W[, 2], lambda = 1))
tt <- seq(0, 1, length.out = 1001)

# Z <- cbind(W)
# plot(Z)

# uncorrected estimators
Acfg.test.m <-
  copula::An.biv(Z, tt, estimator = "CFG", ties.method = "max")
Apick.test.m <-
  copula::An.biv(Z, tt, estimator = "Pickands", ties.method = "max")
Acfg.test.a <-
  copula::An.biv(Z, tt, estimator = "CFG", ties.method = "average")
Apick.test.a <-
  copula::An.biv(Z, tt, estimator = "Pickands", ties.method = "average")
Acfg.test.r <-
  copula::An.biv(Z, tt, estimator = "CFG", ties.method = "random")
Apick.test.r <-
  copula::An.biv(Z, tt, estimator = "Pickands", ties.method = "random")

# true A
A.test <- A(gumbelCopula(3), tt)

# A plot
tmp.test <- eA(x = Z, status = rep(1, n))
# B-splines estimator
m <- fitAspline(Z,
                20,
                tt = tt,
                norder = 4,
                l = 70
)

# Plot for maximum ranks
plot(
  NULL,
  pch = 20,
  ylim = c(0.5, 1),
  xlim = c(0, 1),
  xlab = "t",
  ylab = "A(t)"
)
points(tmp.test[, 1],
       tmp.test[, 2],
       col = "lightgray",
       pch = 20
)
lines(c(0, 0.5, 1),
      c(1, 0.5, 1),
      col = "lightgray",
      lty = 2
)
lines(c(0, 1), c(1, 1), col = "lightgray", lty = 2)

lines(tt, A.test, col = "black")

# CFG estimator
# uncorrected
lines(tt, pmin(1, pmax(Acfg.test.m, tt, 1 - tt)), type = "l", col = "red")
# Deheuvels endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "cfg",
  ties.method = "max",
  rev = TRUE,
  lty = 2,
  col = "red"
)
# Pickands estimator
# uncorrected
lines(tt, pmin(1, pmax(Apick.test.m, tt, 1 - tt)), type = "l", col = "blue")
# Deheuvels endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "pickands",
  ties.method = "max",
  rev = TRUE,
  lty = 2,
  col = "blue",
  madj = 1
)
# Hall & Tajvidi endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "pickands",
  ties.method = "max",
  rev = TRUE,
  lty = 3,
  col = "blue",
  madj = 2
)
# B splines estimator
lines(tt, m$con$A, col = "darkgreen")

# Plot for random ranks
plot(
  NULL,
  pch = 20,
  ylim = c(0.5, 1),
  xlim = c(0, 1),
  xlab = "t",
  ylab = "A(t)"
)
points(tmp.test[, 1],
       tmp.test[, 2],
       col = "lightgray",
       pch = 20
)
lines(c(0, 0.5, 1),
      c(1, 0.5, 1),
      col = "lightgray",
      lty = 2
)
lines(c(0, 1), c(1, 1), col = "lightgray", lty = 2)

lines(tt, A.test, col = "black")

# CFG estimator
# uncorrected
lines(tt, pmin(1, pmax(Acfg.test.r, tt, 1 - tt)), type = "l", col = "red")
# Deheuvels endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "cfg",
  ties.method = "random",
  rev = TRUE,
  lty = 2,
  col = "red"
)
# Pickands estimator
# uncorrected
lines(tt, pmin(1, pmax(Apick.test.r, tt, 1 - tt)), type = "l", col = "blue")
# Deheuvels endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "pickands",
  ties.method = "random",
  rev = TRUE,
  lty = 2,
  col = "blue",
  madj = 1
)
# Hall & Tajvidi endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "pickands",
  ties.method = "random",
  rev = TRUE,
  lty = 3,
  col = "blue",
  madj = 2
)
# B splines estimator
lines(tt, m$con$A, col = "darkgreen")

# Plot for average ranks
plot(
  NULL,
  pch = 20,
  ylim = c(0.5, 1),
  xlim = c(0, 1),
  xlab = "t",
  ylab = "A(t)"
)
points(tmp.test[, 1],
       tmp.test[, 2],
       col = "lightgray",
       pch = 20
)
lines(c(0, 0.5, 1),
      c(1, 0.5, 1),
      col = "lightgray",
      lty = 2
)
lines(c(0, 1), c(1, 1), col = "lightgray", lty = 2)

lines(tt, A.test, col = "black")

# CFG estimator
# uncorrected
lines(tt, pmin(1, pmax(Acfg.test.a, tt, 1 - tt)), type = "l", col = "red")
# Deheuvels endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "cfg",
  ties.method = "average",
  rev = TRUE,
  lty = 2,
  col = "red"
)
# Pickands estimator
# uncorrected
lines(tt, pmin(1, pmax(Apick.test.a, tt, 1 - tt)), type = "l", col = "blue")
# Deheuvels endpoint correction + applying pmin(1, pmax(A, x, 1 - x)) to make sure it's in the triangle
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "pickands",
  ties.method = "average",
  rev = TRUE,
  lty = 2,
  col = "blue",
  madj = 1
)
abvnonpar(
  x = tt,
  data = Z,
  epmar = TRUE,
  add = TRUE,
  plot = TRUE,
  method = "pickands",
  ties.method = "average",
  rev = TRUE,
  lty = 3,
  col = "blue",
  madj = 2
)
# B splines estimator
lines(tt, m$con$A, col = "darkgreen")

# Examining the strange behavior of the Pickands and CFG estimators in the Figure.
# checking what happens at 0.5

# true xi's without ties
xis <- pmin(-log(W[, 1]) / (0.5),
            -log(W[, 2]) / 0.5)
summary(xis)
hist(xis, main = "", xlab = "xi(t)")

# true xi's with ties
xis <- pmin(-log(Z[, 1]) / (0.5),
            -log(ppois(Z[, 2], lambda = 1)) / 0.5)
summary(xis)
hist(xis, main = "", xlab = "xi(t)")

# estimated xi's
Z.ps <- copula::pobs(Z, ties.method = "random")
xis <- pmin(-log(Z.ps[, 1]) / (0.5), -log(Z.ps[, 2]) / 0.5)
summary(xis)
hist(xis, main = "", xlab = "xi(t)")
