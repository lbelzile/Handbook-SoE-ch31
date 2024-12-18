setwd(this.path::here())
library(mev)
library(longevity)
library(ggplot2)
library(patchwork)
library(dplyr)
library(copula)
library(survival)
source("03-helper_functions.R")
source("04-functions-dependence-analysis.R")
# ggplot theme for graphics
theme_set(theme_classic())
generatePlots <- FALSE
set.seed(202312)
# Load data
data("loss", package = "copula")
# Transform survival data to account for metadata and rounding
n <- nrow(loss)
data <- loss |>
  dplyr::rename(rcens = censored) |>
  dplyr::mutate(
    rcens = 1L - rcens, # use survival convention
    Loss = loss,
    ALAE = alae,
    # Scale to K USD for the optimization
    loss = loss / 1e3,
    alae = alae / 1e3,
    limit = limit / 1e3,
    # Create an indicator of rounded values ('clustering')
    intervalcens = ifelse(Loss %% 500 == 0, 0L, 1L),
    time = ifelse(intervalcens == 0L, loss - 0.25, loss),
    time2 = dplyr::case_when(rcens == 0L ~ Inf,
      intervalcens == 0L ~ loss + 0.25,
      .default = loss
    ),
    event = dplyr::case_when(rcens == 0L ~ 0L,
      time == time2 ~ 1L,
      .default = 3L
    ),
    symbol = factor(event, labels = c("right censored", "observed", "rounded"))
  )


df_plot <- data |> dplyr::arrange(by = rcens, decreasing = TRUE)
# Plot bivariate series
g1 <- ggplot() +
  geom_point(
    data = data |> dplyr::filter(event == 1),
    aes(x = Loss, y = ALAE),
    alpha = 0.1,
    # size = 0.5,
    shape = 4
  ) +
  # Right censored
  geom_point(
    data = data |> dplyr::filter(event == 0),
    aes(x = Loss, y = ALAE),
    shape = 1
  ) +
  geom_point(
    data = data |> dplyr::filter(event == 3),
    aes(x = Loss, y = ALAE),
    alpha = 0.25,
    shape = 3
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1e2, 5e2, 1e3, 1e4, 1e5, 5e5),
    labels = c("100", "500", "1K", "10K", "100K", "0.5M")
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(1e2, 5e2, 1e3, 1e4, 1e5, 1e6),
    labels = c("100", "500", "1K", "10K", "100K", "1M")
  ) +
  labs(
    x = "indemnity payment (in USD)",
    y = "allocated loss adjustment expense (in USD)"
  ) +
  theme(legend.position = "bottom")

# Fit the generalized Pareto model accounting for right-censoring
gp_fit_loss <- with(
  data,
  fit_elife(
    time = loss,
    event = rcens,
    type = "right",
    thresh = 1e2,
    family = "gp",
    export = TRUE
  )
)
# Summary of fit of marginal generalized Pareto
summary(gp_fit_loss)
# Quantile-quantile plot on exponential margins
autoplot(gp_fit_loss, which.plot = "exp")
# Fit a generalized Pareto model for the ALAE series
gp_fit_alae <- mev::fit.gpd(
  xdat = data$alae,
  threshold = 20
) # 20K is 0.85666 quantile
# ecdf(lossalae$alae)(5e4)
# Check goodness-of-fit with Q-Q plot
plot(gp_fit_alae)

# Same, but treating any rounded value as interval censored within 250
gp_fit_loss_ic <- fit_elife(
  arguments = data,
  type = "interval2",
  thresh = 1e2,
  family = "gp",
  export = TRUE
)
summary(gp_fit_loss_ic)
ecdf_loss <- np_elife(
  time = data$time,
  time2 = data$time2,
  event = data$event,
  type = "interval2",
  family = "gp",
  export = TRUE
)
ecdf_alae <- ecdf(data$alae)
lambdau <- 1 - c(ecdf_loss$cdf(100), ecdf_alae(20))
# Compute a marginal distribution using semiparametric model (empirical - GP)
marg_dist <- function(x, ecdf, gpmodel) {
  u <- gpmodel$thresh
  qlev <- ecdf(u)
  ifelse(x <= u,
    ecdf(x),
    qlev + (1 - qlev) *
      mev::pgp(
        q = x,
        loc = u,
        scale = coef(gpmodel)["scale"],
        shape = coef(gpmodel)["shape"]
      )
  )
}
# Obtain uniform positions
unif_loss <-
  cbind(
    marg_dist(
      x = data$time,
      ecdf = ecdf_loss$cdf,
      gpmodel = gp_fit_loss
    ),
    marg_dist(
      x = data$time2,
      ecdf = ecdf_loss$cdf,
      gpmodel = gp_fit_loss
    )
  )
unif_alae <- mev::spunif(x = data$alae, thresh = 20)
# Check agreement with uniform (ignoring censoring)
# plot(x = ppoints(nrow(data)), y = unif_loss[,1])


# Create a plot to show the exceedances in more detail
umin <- 0
umax <- qexp(0.9999)
ue <- qexp(1 - lambdau)
df_uplot <- data.frame(
  loss1 = qexp(unif_loss[, 1]),
  loss2 = qexp(unif_loss[, 2]),
  alae = qexp(unif_alae),
  event2 = ifelse(unif_alae < 0.7, 2, 1),
  event1 = ifelse(unif_loss[, 2] < 0.7, 2, data$event)
) |>
  dplyr::filter(loss1 > ue[1] | alae > ue[2])
g2 <- ggplot() +
  geom_rect(
    data = data.frame(
      xmin = umin,
      xmax = umax,
      ymin = umin,
      ymax = umax
    ),
    aes(
      xmin = xmin,
      ymin = ymin,
      ymax = ymax,
      xmax = xmax
    ),
    fill = "grey90"
  ) +
  geom_segment(
    data = data.frame(
      x = umin,
      xend = umax,
      y = qexp(0.7),
      yend = qexp(0.7)
    ),
    aes(
      x = x,
      y = y,
      yend = yend,
      xend = xend
    ),
    linetype = "dashed"
  ) +
  geom_segment(
    data = data.frame(
      y = umin,
      yend = umax,
      x = qexp(0.7),
      xend = qexp(0.7)
    ),
    aes(
      x = x,
      y = y,
      yend = yend,
      xend = xend
    ),
    linetype = "dashed"
  ) +
  geom_rect(
    data = data.frame(
      xmin = umin,
      xmax = ue[1],
      ymin = umin,
      ymax = ue[2]
    ),
    aes(
      xmin = xmin,
      ymin = ymin,
      ymax = ymax,
      xmax = xmax
    ),
    fill = "white"
  ) +
  geom_point(
    df_uplot |> dplyr::filter(event1 == 1, event2 == 1),
    mapping = aes(x = loss1, y = alae),
    shape = 4
  ) +
  # ALAE left-censored
  geom_segment(
    data = df_uplot |> dplyr::filter(event1 == 1, event2 == 2),
    aes(
      y = 0,
      yend = qexp(0.7),
      x = loss1,
      xend = loss1
    ),
    alpha = 0.75,
    color = "grey"
  ) +
  geom_point(
    data = df_uplot |> dplyr::filter(event1 == 1, event2 == 2),
    aes(y = alae, x = loss1),
    alpha = 0.75,
    color = "grey"
  ) +
  # LOSS left-censored
  geom_segment(
    data = df_uplot |> dplyr::filter(event1 == 2, event2 == 1),
    aes(
      y = alae,
      yend = alae,
      x = 0,
      xend = qexp(0.7)
    ),
    alpha = 0.75,
    color = "grey"
  ) +
  geom_point(
    data = df_uplot |> dplyr::filter(event1 == 2, event2 == 1),
    aes(y = alae, x = loss1),
    alpha = 0.75,
    color = "grey"
  ) +
  # Right censoring for loss
  geom_segment(
    data = df_uplot |> dplyr::filter(event1 == 0, event2 == 1),
    aes(
      y = alae,
      yend = alae,
      x = loss1,
      xend = umax
    )
  ) +
  geom_point(
    data = df_uplot |> dplyr::filter(event1 == 0, event2 == 1),
    aes(y = alae, x = loss1),
    shape = 5
  ) +
  # Interval censoring for loss
  geom_segment(
    data = df_uplot |> dplyr::filter(event1 == 3, event2 == 1),
    aes(
      y = alae,
      yend = alae,
      x = loss1,
      xend = loss2
    ),
    linewidth = 2
  ) +
  geom_point(
    data = df_uplot |> dplyr::filter(event1 == 3, event2 == 1),
    aes(y = alae, x = loss1),
    shape = 3
  ) +
  scale_x_continuous(
    limits = c(umin, umax),
    expand = expansion(mult = c(0, 0.05)),
    oob = scales::squish
  ) +
  scale_y_continuous(
    limits = c(umin, umax),
    expand = expansion(mult = c(0, 0.05)),
    oob = scales::squish
  ) +
  labs(
    x = "loss (exponential scale)",
    y = "ALAE (exponential scale)"
  )

## Likelihood ratio test for equality of shape.
loglik_full <- as.numeric(logLik(gp_fit_loss_ic)) +
  as.numeric(logLik(gp_fit_alae))
# Create a wrapper for the log likelihood for jt estimation
nll_pooled <- function(param) {
  as.numeric(
    mev:::gpd.ll.optim(
      par = c(log(param[2]), param[3]),
      dat = gp_fit_alae$exceedances
    ) +
      longevity::nll_elife(
        par = c(param[1], param[3]),
        time = data$time,
        time2 = data$time2,
        event = data$event,
        type = "interval2",
        thresh = 1e2,
        family = "gp"
      )
  )
}
# Use fit from individual observations
start <- c(coef(gp_fit_loss_ic)["scale"], coef(gp_fit_alae))
# Check that the starting value is feasible
nll_pooled(param = start)
# Compute the pooled mle
opt_pooled <-
  optim(
    par = start,
    fn = nll_pooled,
    method = "Nelder-Mead",
    hessian = TRUE
  )
# Calculate standard errors
sqrt(diag(solve(opt_pooled$hessian)))
# Compute the p-value for the likelihood ratio test using asymptotic null dist.
pchisq(2 * (opt_pooled$value + loglik_full), df = 1, lower.tail = FALSE)


# Compute tail correlation coefficient
# accounting for survival problems
#
# We want Pr(F_1(Y_1) > u, F_2(Y_2) > u)/(1-u)
# Usual estimator is Pr(min_j F_j(Y_j) > u)/(1-u)
#
# Also could write this as Pr(F_1(Y_1) > u | F_2(Y_2) > u)
# Below, a poor man's attempt at retrieving this quantity
# by first fitting nonparametric MLE to obtain F_1^{-1}(u),
# then estimating the survival probability above that level for LOSS (Y_1)
# conditional on exceedances of ALAE (Y_2)
if (generatePlots) {
  unif_alae_np <- rank(data$alae) / (n + 1)
  qlev <- seq(0.8, 0.98, by = 0.01)
  npsurv_loss <- np_elife(
    type = "interval2",
    thresh = 0,
    arguments = data
  )
  quants <- quantile(
    x = npsurv_loss$cdf,
    prob = qlev
  )
  chi <- numeric(length(qlev))
  for (i in seq_along(qlev)) {
    chi[i] <- (1 - np_elife(
      type = "interval2",
      thresh = 0,
      arguments = data[unif_alae_np > qlev[i], ]
    )$cdf(quants[i]))
  }
  # An alternative strategy would be to map data to uniform using
  # the nonparametric maximum likelihood estimator
  # then to censor the variable accordingly and
  struct_unif <- pmin(npsurv_loss$cdf(data$time), unif_alae_np)
  struct_unif2 <- pmin(npsurv_loss$cdf(data$time2), unif_alae_np)
  np_struct_unif <- np_elife(
    time = struct_unif,
    time2 = struct_unif2,
    event = case_when(struct_unif != struct_unif2 ~ 3L,
      .default = 1L
    ),
    type = "interval2"
  )
  pexc <- 1 - np_struct_unif$cdf(qlev)
  g4 <- ggplot() +
    geom_pointrange(
      data = data.frame(
        qlev = qlev,
        chi = pexc / (1 - qlev)
      ),
      mapping = aes(
        x = qlev,
        y = chi,
        ymin = chi + qnorm(0.025) * sqrt(pexc * (1 - pexc)) / (1 - qlev) /
          sqrt(n),
        ymax = chi + qnorm(0.975) * sqrt(pexc * (1 - pexc)) / (1 - qlev) /
          sqrt(n)
      ),
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    labs(
      x = "quantile level",
      y = "",
      subtitle = expression("Tail correlation" ~ chi)
    )

  # Compute eta through MLE
  eta <- matrix(nrow = length(qlev), ncol = 3)
  for (i in seq_along(qlev)) {
    fit_exp <- fit_elife(
      time = qexp(struct_unif),
      time2 = qexp(struct_unif2),
      event = case_when(struct_unif != struct_unif2 ~ 3L,
        .default = 1L
      ),
      type = "interval2",
      thresh = qexp(qlev[i]),
      family = "exp"
    )
    eta[i, 1] <- coef(fit_exp)
    eta[i, -1] <- confint(fit_exp)
  }
  g3 <- ggplot(
    data = data.frame(
      qlev = qlev,
      eta = eta[, 1],
      lower = eta[, 2],
      upper = eta[, 3]
    ),
    mapping = aes(
      x = qlev,
      y = eta,
      ymin = lower,
      ymax = upper
    )
  ) +
    geom_pointrange() +
    scale_y_continuous(
      limits = c(0, 1),
      oob = scales::squish,
      expand = c(0, 0),
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    labs(
      x = "quantile level",
      y = "",
      subtitle = expression("Tail dependence" ~ eta)
    )
  g3 + g4
  ggsave(
    filename = "../LaTeX/figures/lossalae_fig2.pdf",
    width = 8,
    height = 4
  )
}

# Generate Q-Q plot pooling two series to save space
if (generatePlots) {
  loss_res1 <- mev::pgp(
    q = gp_fit_loss_ic$time,
    scale = coef(gp_fit_loss_ic)["scale"],
    shape = coef(gp_fit_loss_ic)["shape"]
  )
  loss_res2 <- mev::pgp(
    q = gp_fit_loss_ic$time2,
    scale = coef(gp_fit_loss_ic)["scale"],
    shape = coef(gp_fit_loss_ic)["shape"]
  )
  alae_res <- mev::pgp(
    q = gp_fit_alae$exceedances,
    scale = coef(gp_fit_alae)["scale"],
    shape = coef(gp_fit_alae)["shape"]
  )
  xpos_loss <- qexp(
    np_elife(
      time = gp_fit_loss_ic$time,
      time2 = gp_fit_loss_ic$time2,
      event = gp_fit_loss_ic$event,
      type = "interval2"
    )$cdf((
      gp_fit_loss_ic$time + gp_fit_loss_ic$time2
    ) / 2)
  )
  # Use average positions ('clustered' data positions are the ones reported)
  xpos_alae <-
    qexp(rank(gp_fit_alae$exceedances) / (gp_fit_alae$nat + 1))
  nexc_pooled <- length(xpos_loss) + length(xpos_alae)
  set.seed(2024)
  B <- 1000L
  mmax <- 7
  # Compute tolerance bands via Monte Carlo (could use Beta order stats, this is faster)
  tolbands <- t(apply(
    apply(
      matrix(
        rexp(B * 1000),
        nrow = B
      ), 1,
      sort
    ), 1,
    quantile,
    probs = c(0.05, 0.95)
  ))
  g6 <- ggplot() +
    geom_ribbon(
      data = data.frame(
        x = qexp(ppoints(1000)),
        ymin = tolbands[, 1],
        ymax = tolbands[, 2]
      ),
      aes(x = x, ymin = ymin, ymax = ymax),
      fill = "gray",
      alpha = 0.25
    ) +
    geom_point(
      data = data.frame(x = xpos_alae, y = qexp(alae_res)),
      aes(x = x, y = y)
    ) +
    geom_point(
      data = data.frame(x = xpos_loss, y = qexp(loss_res1)) |>
        dplyr::filter(gp_fit_loss_ic$event %in% c(1, 3)),
      aes(x = x, y = y)
    ) +
    geom_segment(
      data = data.frame(
        x = xpos_loss,
        y = qexp(loss_res1),
        yend = qexp(loss_res2)
      ) |>
        dplyr::filter(gp_fit_loss_ic$event == 3),
      aes(
        x = x,
        y = y,
        xend = x,
        yend = yend
      ),
      linewidth = 3,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0, mmax),
      expand = expansion(mult = c(0, 0.05)),
      oob = scales::squish
    ) +
    scale_y_continuous(
      limits = c(0, mmax),
      expand = expansion(mult = c(0, 0.05)),
      oob = scales::squish
    ) +
    labs(y = "standardized observed quantiles", x = "exponential quantiles") +
    theme_classic()
}


if (generatePlots) {
  g1 + g6
  ggsave("../LaTeX/figures/lossalae_fig1.pdf",
         width = 8,
         height = 4
  )
}

###########################################################
####  Multivariate generalized Pareto logistic model   ####
###########################################################
# Starting parameters:  divide scale by 100
# (along with observations, to make parameters comparable)
start <- par <- c(opt_pooled$par / c(100, 100, 1), alpha = 0.7)
## Fit a multivariate logistic model
# thresholds in K dollars, respectively ~90% and ~95% marginal quantiles
# Fit the model with marginal left-censoring, but using the MEV copula
# but the model uses information below, so needs all observations...
opt_fbvpot <- evd::fbvpot(
  cbind(data$loss, data$alae) / 100,
  threshold = c(100, 20) / 100,
  model = "log",
  likelihood = "censored",
  start = list(
    scale1 = start[1],
    scale2 = start[2],
    shape1 = start[3],
    dep = 0.5
  ),
  cshape = TRUE
)


######################################################
# USING MEV LIKELIHOOD, ACCOUNTING FOR RIGHT-CENSORING
# Check for marginal fit via threshold stability plots
longevity::tstab(
  time = data$time,
  time2 = data$time2,
  event = data$event,
  thresh = quantile(ecdf_loss$cdf, seq(0.6, 0.98, by = 0.02)),
  family = "gp",
  type = "interval2"
)
plot(
  longevity::nc_test(
    time = data$time,
    time2 = data$time2,
    event = data$event,
    thresh = quantile(ecdf_loss$cdf,
                      probs = seq(0.7, 0.98, by = 0.02)),
    type = "interval2"
  )
)
mev::tstab.gpd(
  data$alae,
  thresh = quantile(data$alae,
                    probs = seq(0.6, 0.98, by = 0.02)))
plot(longevity::nc_test(
  time = data$alae,
  thresh = quantile(data$alae,
                    probs = seq(0.6, 0.98, by = 0.02))))
# It appears that the generalized Pareto model could fit
# at lower level for ALAE than for losses
thresh <- c(100, 20)
plot(fit.gpd(data$alae, threshold = 20))
ecdf_loss$cdf(100)
ecdf_alae(20)


# exceeds <- which(data$loss > thresh[1] | data$alae > thresh[2])
# Ntot <- nrow(data)
# # Percentage of exceedances
# nexc <- length(exceeds)
# (pat <- nexc / Ntot)
# Marginal threshold for inferential left-censoring
# given as a probability
margpcens <- rep(0.5, 2L)
# Transform data into an array
adat <- array(NA, dim = c(nrow(data), 2, 2))
adat[,1,] <- cbind(data$time, data$time2)
adat[,2,] <- cbind(data$alae, data$alae)
# Create vector of bivariate indicators
biv_event <- cbind(as.integer(data$event), rep(1L, nrow(data)))
nll_pot_logist <- function(par,
                        margpcens = c(0.5,0.5),
                        margthresh = c(100, 20),
                        depthresh = 0.8,
                        cshape = TRUE,
                        profile = NULL,
                        fpar = NA
                        ){
  if(!is.null(profile)){
    if(profile == "shape"){
      par <- c(par[1:2], fpar, par[-(1:2)])
    } else if(profile == "alpha"){
      par <- c(par, fpar)
    }
  }
  lpar <- ifelse(isTRUE(cshape), 4L, 5L)
  stopifnot(length(par) == lpar)
  scale <- par[1:2]
  if(isTRUE(cshape)){
   shape <- rep(par[3], 2)
   alpha <- par[4]
  } else{
    shape <- par[3:4]
    alpha <- par[5]
  }

  # Some sanity constraints for the parameters
  if(isTRUE(any(scale <= 0.05, shape <= -1, shape > 1, alpha < 0))){
    return(1e10)
  }
  -bvpot_log(dat = adat,
             event = biv_event,
             depthresh = depthresh,
             margthresh = margthresh,
             margpcens = margpcens,
             scale = scale,
             shape = shape,
             par = alpha,
             mdist = "mev",
             ecdf1 = ecdf_loss$cdf,
             ecdf2 = ecdf_alae)
}
# Check that starting values are feasible
start <- opt_fbvpot$estimate * c(100, 100, 1, 1)
cshape <- TRUE
profile <- function(par = c("alpha","shape"),
                    level = 0.95,
                    thresh,
                    margpcens,
                    start,
                    ...){
  cshape <- TRUE
  par <- match.arg(par)
  opt <- optim(
    fn = nll_pot_logist,
    par = start,
    depthresh = thresh,
    margpcens = margpcens,
    cshape = cshape,
    method = "Nelder",
    control = list(parscale = c(rep(100, 2), rep(1, 2 + !cshape)),
                   reltol = 1e-8,
                   abstol = 1e-10,
                   maxit = 2000),
    hessian = TRUE)
  mle <- opt$par
  if(mle[4] > 1){
    mle[4] <- 1/mle[4]
  }
  se <- 0.08
  fparind <- switch(par, alpha = 4L, shape = 3L) + !cshape
  fpar_grid <- mle[fparind] + seq(-2.5*se, 2.5*se, length.out = 31)
  fpar_grid <- fpar_grid[fpar_grid > 0 & fpar_grid < 1]
  nll <- numeric(length = length(fpar_grid))

  for(i in seq_along(fpar_grid)){
    opti <- try(nlminb(
    objective = nll_pot_logist,
    start = mle[-fparind],
    depthresh = thresh,
    profile = par,
    fpar = fpar_grid[i],
    margpcens = margpcens,
    lower = c(0.01, 0.01,-0.99,0.01),
    upper = c(Inf, Inf, 1.5, 0.99)))
    if(!inherits(opti, "try-error")){
      nll[i] <- opti$objective
    }
    # control = list(parscale = c(rep(100, 2), rep(1, 1 + !cshape)),
    #                reltol = 1e-8, abstol = 1e-10))$value
  }
  prof <- list(pll = -nll, psi.max = mle[fparind], maxpll = -opt$value, psi = fpar_grid)
  mev:::confint.eprof(object = prof, parm = "profile", level = level,method = "smooth.spline")
}
bvpot_opt <- optim(
  fn = nll_pot_logist,
  par = start,
  method = "BFGS",
  control = list(parscale = c(rep(100, 2), rep(1, 2 + !cshape)),
                 reltol = 1e-8, abstol = 1e-10),
  hessian = TRUE)
numDeriv::grad(nll_pot_logist,
               x = bvpot_opt$par)
# Did we converge to an optimum?
min(eigen(bvpot_opt$hessian, only.values = TRUE)$values) > 0
# Parameter estimates and standard errors
mle <- bvpot_opt$par
se <- sqrt(diag(solve(bvpot_opt$hessian)))
(profile("shape",
        thresh = 0.8,
        margpcens = c(0.5,0.5),
        start = mle))

#### WARNING - COMPUTATIONALLY INTENSIVE
# Sensitivity analysis
tlev <- seq(0, 0.9, by = 0.1)
mlev <- c(0, 0.25, 0.5, 0.75)
sensmat <- matrix(nrow = 0, ncol = 5)
ind <- 1
for(i in seq_along(tlev)){
  for(j in seq_along(mlev)){
    if(mlev[j] <= tlev[i]){
      prof_dep <- try(profile("alpha", thresh = tlev[i], margpcens = rep(mlev[j], 2), start = mle))
      if(inherits(prof_dep, "try-error")){
        prof_dep <- rep(NA, 3)
      }
      # Paper uses theta = 1/alpha, but easier to optimize over alpha
      sensmat <- rbind(sensmat, c(tlev[i], mlev[j], 1/prof_dep))
      ind <- ind + 1L
    }
  }
}
# Transform to data frame and cast variables to factors
colnames(sensmat) <- c("depthresh", "margcensthresh", "mle", "upper", "lower")
sens_dat <- as.data.frame(sensmat) |>
  dplyr::mutate(margcensthresh = factor(margcensthresh),
                depthresh = factor(depthresh))
g7 <- ggplot(data = sens_dat,
             mapping = aes(x = depthresh,
                           group = margcensthresh,
                           color = margcensthresh,
                           y = mle,
                           ymin = lower,
                           ymax = upper)) +
  geom_pointrange(position = position_dodge(width = 0.33),
                  size = 0.1) +
  scale_y_continuous(limits = c(1,1.6),
                     expand = c(0,0),
                     breaks = seq(1, 1.6, by = 0.1),
                     labels = c("1","1.1","1.2","1.3","1.4","1.5","1.6")) +
  scale_color_grey() +
  labs(y = expression(theta),
       x = "dependence threshold (prob. scale)",
       color = "marg. censoring level") +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey90",
                                          linewidth = 0.1))


#### Dependence analysis using copulas

# CAREFUL: code is restricted to the case discussed in the book,
# i.e. only one variable is censored and has ties.

# data are as calculated on line 22
# cop_data contains the variables used in Section 1.3.3 :
# loss (1st column), ALAE (2nd column),
# and censoring status of loss (3rd column)

set.seed(202312)
cop_data <- data |>
  dplyr::select(c("Loss", "ALAE", "rcens"))

# code for Figure 1.6 (left panel)

# calculating the A-plot points

tmp <- eA(x = cop_data[, c(1, 2)], status = cop_data[, 3])

# Estimating Pickands dependence function using the Pickands and CFG estimators, ignoring ties and censoring

tt <- seq(0, 1, length.out = 1001)
Acfg <- copula::An.biv(cop_data[, c(1, 2)], tt, estimator = "CFG")
Apick <-
  copula::An.biv(cop_data[, c(1, 2)], tt, estimator = "Pickands")

# Calculating the spline estimator from Cormier et al (2014) and Bücher et al (2022)

# selecting the penalty with LOO-cross validation
# WARNING: computationally intensive
lambdas <- seq(5L, 75L, by = 5L)
cvs <- fitAsplineCV(
  x = cop_data[, c(1, 2)],
  lambdas,
  l = 40,
  norder = 4,
  status = cop_data[, 3]
)
# The grid of lambda is too small, and increasing the penalty
# keeps reducing the error. However, the differences between
# the fitted A curves are indistinguishable to the naked eye.
plot(lambdas,
  cvs,
  type = "b",
  xlab = "lambda",
  ylab = "CV"
)

# calculating the B-spline estimator

m4 <- fitAspline(
  cop_data[, c(1, 2)],
  lambda = 15,
  tt = tt,
  norder = 4,
  l = 40,
  status = cop_data[, 3]
)

# spline of order 3, not used in the book
# m3 <- fitAspline(
#  cop_data[,c(1,2)],
#  15,
#  tt = tt,
#  norder = 3,
#  l = 40,
#  status = cop_data[,3])

# plotting
if (generatePlots) {
  g5 <- ggplot(
    data = data.frame(
      t = tmp[, 1],
      A = tmp[, 2]
    ),
    mapping = aes(x = t, y = A)
  ) +
    geom_polygon(
      data = data.frame(
        x = c(0, 0.5, 1, 0),
        y = c(1, 0.5, 1, 0)
      ),
      mapping = aes(x = x, y = y),
      linetype = "dashed",
      fill = "white",
      color = "gray"
    ) +
    geom_point(alpha = 0.25, col = "gray80") +
    geom_line(data = data.frame(t = tt, A = Acfg), linetype = "dashed") +
    geom_line(data = data.frame(t = tt, A = Apick), linetype = "dotted") +
    geom_line(data = data.frame(t = tt, A = m4$con$A)) +
    scale_y_continuous(
      limits = c(0.5, 1),
      breaks = c(0.5, 0.75, 1),
      labels = c("0.5", "0.75", "1"),
      expand = expansion(mult = c(0.05))
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      expand = expansion(mult = c(0.05))
    ) +
    labs(
      y = "Pickands dependence function, A(w)",
      x = "angle w"
    )

  g5 + g2
  ggsave(
    file = "../LaTeX/figures/lossalae_fig3.pdf",
    width = 8,
    height = 4
  )
}

## Applying the tests of extremeness

evTestK(as.matrix(cop_data[, c(1, 2)]),
        method = "asymptotic",
        ties = TRUE)

# Test of bivariate extreme-value dependence based on Kendall's
# distribution with argument 'method' set to "asymptotic"
# statistic = -0.59774, p-value = 0.55

evTestA(as.matrix(cop_data[, c(1, 2)]),
        derivatives = "Cn",
        ties.method = "max")
# 	Test of bivariate extreme-value dependence based on the CFG
# 	estimator with argument 'derivatives' set to 'Cn'
# statistic = 0.019039, p-value = 0.5579

evTestA(as.matrix(cop_data[, c(1, 2)]),
        derivatives = "Cn",
        ties.method = "random")

# Test of bivariate extreme-value dependence based on the CFG
# estimator with argument 'derivatives' set to 'Cn'
# statistic = 0.017265, p-value = 0.461

## Estimating the parameter of the Gumbel--Hougaard copula

X <- cop_data[, 1] # losses
Y <- cop_data[, 2] # alae
status <- cop_data[, 3] # censoring status

# calculating the pseudos (Kaplan-Meier for losses, ECDF for ALAE)
pseudos <- GetPseudosEst(X, Y, status)

# initial parameter value calculated from the likelihood ignoring ties and censoring
theta.init <-
  summary(fitCopula(gumbelCopula(), pobs(cbind(
    pseudos$U, pseudos$V
  )), method = "mpl"))$coefficients[1]

# Likelihood with the Ery tie adjustment
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
  upper = 2
)

# $par
# [1] 1.438439

# $value
# [1] 18258.63


# Likelihood with the Li et al tie adjustment
optim(
  theta.init,
  fn = negloglik.cens.Li,
  U = pseudos$U,
  U.left = pseudos$U.left,
  V = pseudos$V,
  V.left = pseudos$V.left,
  status = status,
  tie = pseudos$tie,
  method = "Brent",
  lower = 1,
  upper = 2
)

# $par
# [1] 1.438491
#
# $value
# [1] 4511.691

## I observed that sometimes Brent does really poorly, even though the likelihood looks beautiful when plotted. Nelder-Mead may be more stable.

# What follows is not included in the book: to check, I extracted the uncensored observations
# then used the code from Nasri & Rémillard (library(CopulaInference)) with the Li et al. tie adjustment
# I also estimated theta by ignoring ties and using the code from the copula package with the three tie breaking methods.

# extracting uncensored observations
library(CopulaInference)
cop_data.u <- cop_data[cop_data[, 3] == 1, ]
out.u <-  EstBiCop(data = cop_data.u[, c(1, 2)], family = "gumbel")
# From Nasri & Remillard
summary(out.u)
# out.u$par
# [1] 1.424866

fitCopula(gumbelCopula(), pobs(cop_data.u[, c(1, 2)], ties.method = "max"))
# Call: fitCopula(gumbelCopula(), data = pobs(cop_data.u))
# Fit based on "maximum pseudo-likelihood" and 1466 2-dimensional observations.
# Copula: gumbelCopula
# alpha
# 1.428
# The maximized loglikelihood is 191.4
# Optimization converged

fitCopula(gumbelCopula(), pobs(cop_data.u[, c(1, 2)], ties.method = "average"))
# Call: fitCopula(gumbelCopula(), data = pobs(cop_data.u))
# Fit based on "maximum pseudo-likelihood" and 1466 2-dimensional observations.
# Copula: gumbelCopula
# alpha
# 1.425
# The maximized loglikelihood is 190.9
# Optimization converged

fitCopula(gumbelCopula(), pobs(cop_data.u[, c(1, 2)], ties.method = "random"))
# Call: fitCopula(gumbelCopula(), data = pobs(cop_data.u))
# Fit based on "maximum pseudo-likelihood" and 1466 2-dimensional observations.
# Copula: gumbelCopula
# alpha
# 1.425
# The maximized loglikelihood is 190.7
# Optimization converged

### Running the goodness of fit procedures

# Calculating the ECDF of ALAE and the Kaplan--Meier estimator of loss

km <- survfit(Surv(cop_data[, 1], cop_data[, 3]) ~ 1)
x <- as.numeric(km$time)
Fx <-  as.numeric(1 - km$surv)
# for the loss data, the largest variable is uncensored
# so the KM estimator is 1 at that observation

y <- sort(unique(cop_data[, 2]))
G <- edfun::edfun(cop_data[, 2])
Fy <- G$pfun(y)

# Calculating the policy limit

data(loss)
limit <- loss$limit
# policylimit <- (3 * 10 ^ 6) * (cop_data[, 3] == 1) + cop_data[, 1] * (cop_data[, 3] == 0)
policylimit <-
  (10 * 10^6) * (limit == -99) + limit * (limit > -99)
inds <- order(cop_data[, 2])[!duplicated(sort(cop_data[, 2]))]
p <- policylimit[inds]
# p <- policylimit[order(rank(cop_data[,2]))]

# The goodness-of-fit procedure at the end of Section 1.3.3.2

# the estimate with the Ery tie correction
theta <- 1.438439
# theta <- 1.438491 #Li et al

# calculating the observed value of the test statistic
myRn <- Rn(cop_data[, c(1, 2)], status = cop_data[, 3])

# doing the bootstrap
# THIS IS SUPER SLOW, the loop can be eliminated with library(simsalapar)

N <- 1000
Rn.boot <- numeric(N)
for (i in 1:N) {
  if (i / 100 == round(i / 100)) {
    print(i)
  }
  bootdata <-
    sample.boot(
      1500,
      theta = theta,
      seed = i,
      x = x,
      Fx = Fx,
      y = y,
      Fy = Fy,
      p = p
    )
  Rn.boot[i] <- Rn(bootdata[, c(1, 2)], status = bootdata[, 3])
}

# p-value

mean(Rn.boot > myRn)

# 0.177

# There were a handful of really large values in Rn.boot (5 out of 1000). I investigated this and it is due to unstable optimization of the likelihood when estimating theta (it's the Brent algorithm). But the likelihood curve looks beautiful, so I'm not sure why optim does some nonsence.


# The goodness-of-fit procedure at the end of Section 1.3.3.1

# calculating the observed value of the test statistic
myRn.n <- Rn.n(cop_data[, c(1, 2)],
  status = cop_data[, 3],
  lambda = 15,
  l = 40
)

# doing the bootstrap
# THIS IS AGAIN SUPER SLOW; this could be run using parallel computing

N <- 1000
Rn.n.boot <- numeric(N)
for (i in 1:N) {
  if (i / 100 == round(i / 100)) {
    print(i)
  }
  bootdata <-
    sample.boot.n(
      cop_data[, c(1, 2)],
      15,
      status = cop_data[, 3],
      norder = 4,
      seed = i,
      x = x,
      Fx = Fx,
      y = y,
      Fy = Fy,
      p = p,
      l = 40
    )
  Rn.n.boot[i] <-
    Rn.n(bootdata[, c(1, 2)],
      status = bootdata[, 3],
      lambda = 15,
      l = 40
    )
}

# p-value

mean(Rn.n.boot > myRn.n)

# 0.251
# Assessing the standard error of theta-hat with the Ery correction

theta <- 1.438439

N <- 1000
theta1.boot <- numeric(N)
for (i in 1:N) {
  if (i / 100 == round(i / 100)) {
    print(i)
  }
  bootdata <-
    sample.boot(
      1500,
      theta = theta,
      seed = i,
      x = x,
      Fx = Fx,
      y = y,
      Fy = Fy,
      p = p
    )
  pseudos <- GetPseudosEst(bootdata[, 1], bootdata[, 2], bootdata[, 3])
  theta1.boot[i] <- suppressWarnings(
    optim(
      theta.init,
      fn = negloglik.cens.Ery,
      U = pseudos$U,
      U.left = pseudos$U.left,
      V = pseudos$V,
      V.left = pseudos$V.left,
      status = bootdata[, 3],
      method = "Brent",
      lower = 1,
      upper = 2
    )$par
  )
}

# Assessing the standard error of theta-hat with the Li et al correction

theta <- 1.438491 # Li et al

N <- 1000
theta2.boot <- numeric(N)
for (i in 1:N) {
  if (i / 100 == round(i / 100)) {
    print(i)
  }
  bootdata <-
    sample.boot(
      1500,
      theta = theta,
      seed = i,
      x = x,
      Fx = Fx,
      y = y,
      Fy = Fy,
      p = p
    )
  pseudos <- GetPseudosEst(bootdata[, 1], bootdata[, 2], bootdata[, 3])
  theta2.boot[i] <- suppressWarnings(
    optim(
      theta.init,
      fn = negloglik.cens.Li,
      U = pseudos$U,
      U.left = pseudos$U.left,
      V = pseudos$V,
      V.left = pseudos$V.left,
      status = bootdata[, 3],
      tie = pseudos$tie,
      method = "Brent",
      lower = 1,
      upper = 2
    )$par
  )
}



#####################################################
#### Simulation from the total repayment amount  ####
#####################################################

# Quantile function for the semiparametric marginal model
qsemipar <- function(x, zetau, ecdf, thresh, gppar) {
  out <- numeric(length = length(x))
  xlow <- x <= 1 - zetau
  out[xlow] <- quantile(ecdf, x[xlow])
  out[!xlow] <-
    mev::qgp(
      p = x[!xlow],
      loc = thresh,
      scale = gppar[1],
      shape = gppar[2]
    )
  return(out)
}

# We simulate from the bivariate logistic copula
nsim <- 1e6
opt <- optim(
  fn = nll_pot_logist,
  margpcens = c(0,0),
  depthresh = 0,
  par = start,
  method = "BFGS",
  control = list(parscale = c(rep(100, 2), rep(1, 2 + !cshape)),
                 reltol = 1e-8, abstol = 1e-10),
  hessian = TRUE)
mle <- opt$par
vcov <- solve(opt$hessian)
# Total value of reinsurance
total_reins <- function(loss, alae, r, z) {
  X <- pmin(loss, z)
  mean(ifelse(X > r, X - r + (X - r) / X * alae, 0))
}
# Compute the expected value of total reinsurer payment
r <- seq(20, 150, by = 10)
nr <- length(r)
# Threshold and quantile level of the latter
thresh <- c(100, 20)
zetau <- 1 - c(ecdf_loss$cdf(100), ecdf_alae(20))
nrep <- 100
tot_reins <- matrix(nrow = nrep, ncol = 2*nr)
for(i in seq_len(nrep)){
  set.seed(i)
  print(i)
  # Simulate parameters from multinormal approximation to account
  # for parameter uncertainty
  pars <- mev::mvrnorm(n = 1, mu = mle, Sigma = vcov)
  # Simulate from logistic copula over entire domain
 unif_samp <- copula::rCopula(
  n = nsim,
  copula::archmCopula(
    family = "gumbel",
    param = 1/pars[4],
    dim = 2
  )
)
# Sample policy caps with replacement (this is a proxy for maximum value)
limits_samp <- sample(data$limit, size = nsim, replace = TRUE)

loss_samp <- qsemipar(
  x = unif_samp[, 1],
  zetau = zetau[1],
  ecdf = ecdf_loss$cdf,
  thresh = thresh[1],
  gppar = pars[c(1,3)]
)
alae_samp <- qsemipar(
  x = unif_samp[, 2],
  zetau = zetau[2],
  ecdf = ecdf(data$ALAE),
  thresh = 20,
  gppar = pars[2:3]
)
for (j in seq_along(r)) {
  # Compute the total composite variable
  tot_reins[i,j] <-
    total_reins(
      loss = loss_samp,
      alae = alae_samp,
      z = limits_samp,
      r = r[j]
    )
  # Monte Carlo std.errors are negligible, so not computed
  # Same, but making the policy cap infinite
  tot_reins[i,j+nr] <-
    total_reins(
      loss = loss_samp,
      alae = alae_samp,
      z = Inf,
      r = r[j]
    )
}
}
g8 <- ggplot(data = data.frame(
      x = c(rep(r, each = nrep), rep(r, each = nrep)),
      y = c(tot_reins),
      cap = factor(rep(c(
        "policy limits", "infinite"
      ), each = length(tot_reins)/2))
    ),
    aes(
      x = x,
      y = y,
      color = cap,
      group = cap
    )
  ) +
  ggdist::stat_pointinterval(point_interval = "median_qi",
                             .width = c(0.5, 0.95),
                             point_size = 0.8) +
  ggplot2::scale_color_grey() +
  labs(
    x = "average insurer retention r (in USD)",
    y = "reinsurer premium T (in USD)"
  ) +
  scale_x_continuous(labels = scales::label_number(suffix = " K", scale = 1)) +
  scale_y_continuous(
    labels = scales::label_number(suffix = " M", scale = 1e-3),
    limits = c(0, 3e3),
    expand = c(0, 0)
  ) + # data already in thousands
  theme(legend.position = "bottom")
g7 + g8
ggsave(
  filename = "../LaTeX/figures/lossalae_fig4.pdf",
  width = 10,
  height = 5
)
