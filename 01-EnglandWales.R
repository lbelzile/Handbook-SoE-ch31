setwd(this.path::here())
# Download latest version of package
library(longevity) # extremes and survival analysis
library(dplyr) # data manipulation
library(ggplot2) # grammar of graphics
library(lubridate) # date manipulation
library(patchwork) # combine ggplot objects
library(scales) # scales for plots (dates, etc.)
library(mgcv) # Smoothing
library(gratia) # Methods for GAMS (plots and simultaneous intervals)
theme_set(theme_classic())
generateFigure <- FALSE # binary trigger for figures
# Load IDL data for England and Wales
load("../Data/englandwales.rda")

# The data already includes two columns with left and right truncation bounds
# Counts of occurrence from pooled ONS and IDL data
englandwales |>
  dplyr::mutate(group = factor(ageyear >= 110,
    labels = c("[105, 110)", "[110+")
  )) |>
  dplyr::group_by(group, gender) |>
  dplyr::summarize(count = dplyr::n())
# Check collection period
englandwales |>
  dplyr::mutate(group = factor(ageyear >= 110,
    labels = c("[105, 110)", "[110+")
  )) |>
  dplyr::group_by(group) |>
  dplyr::summarize(
    begin = min(ddate),
    end = max(ddate)
  )
# Counts per age at death, rounded down
englandwales |>
  dplyr::group_by(ageyear) |>
  dplyr::summarize(count = dplyr::n())
# Create a copy of the database with data in years
ew <- with(englandwales, list(
  time = ndays / 365.25,
  ltrunc = cbind(ltrunc1, ltrunc2) / 365.25,
  rtrunc = cbind(rtrunc1, rtrunc2) / 365.25,
  byear = lubridate::year(bdate),
  gender = gender
))


####################################################
###   Life table approach and smoothed hazard    ###
####################################################
tyr_smooth <- 105
thresh_smooth <- ceiling(tyr_smooth * 365.25)
# Bin data to account for truncation
rdays <- c(thresh_smooth, max(englandwales$ndays))
bins <- seq(
  from = rdays[1],
  to = rdays[2],
  length.out = diff(rdays)
)

# Number at risk (conditional on NOT having failed)
nrisk <- sapply(bins, function(bin) {
  with(
    englandwales,
    sum(ltrunc1 <= bin & bin < rtrunc1 & ndays >= bin) +
      sum(ltrunc2 <= bin & bin < rtrunc2 & ndays >= bin, na.rm = TRUE)
  )
})
nonempty <- which(nrisk > 0)
nrisk <- nrisk[nonempty]
# Number of failures
bin_death <- (englandwales$ndays - rdays[1])
# Scale time to start at threshold
midbin <- seq.int(diff(rdays))[nonempty] / 365.25
nfail <- as.integer(table(factor(bin_death, levels = 1:max(bin_death)))[nonempty])


# Create data frame for Poisson table
pois_dat <- data.frame(
  time = (as.integer(bins) - thresh_smooth),
  nrisk = as.integer(nrisk),
  nfail = as.integer(nfail),
  midbin = midbin
)[nrisk > 0, ]

# Switch perspective and use semiparametric model with smoothing P-spline
smoothhazard105 <- mgcv::gam(
  nfail ~ offset(log(nrisk)) + s(midbin, bs = "ps", k = 10),
  family = poisson,
  data = pois_dat,
  method = "REML"
)
smoothhazard105_GCV <- mgcv::gam(
  nfail ~ offset(log(nrisk)) + s(midbin, bs = "ps", k = 10),
  family = poisson,
  data = pois_dat
)
# The default, generalized cross-validation, appears to be oversmoothing
# Test conclusions are the same

# Simultaneous confidence intervals discussed in this post by Gavin Simpson
# https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
conf_interv105 <- gratia:::confint.gam(
  smoothhazard105,
  transform = exp,
  unconditional = TRUE,
  parm = "s(midbin)",
  type = "simultaneous",
  shift = TRUE
)
ci105_GCV <- gratia:::confint.gam(
  smoothhazard105_GCV,
  transform = exp,
  parm = "s(midbin)",
  type = "simultaneous",
  shift = TRUE
)
# Plot the smooth for REML
plot_gam <- plot(smoothhazard105,
  trans = exp,
  rug = FALSE
)
####################################################
###        Peaks over threshold analysis         ###
####################################################
## Threshold stability plot
tstab <- longevity::tstab(
  arguments = ew,
  thresh = seq(105, 110, by = 0.5),
  family = "gp",
  plot.type = "ggplot",
  plot = FALSE,
  method = "profile"
)

# Pick 108 as threshold
tyr <- 108 # threshold (in years)
thresh <- ceiling(365.25 * tyr) # 38352 days

# Keep only plot with shape parameters and fix label
# g2 <- plot(tstab, plot.type = "ggplot")$g1 +
# labs(x = "threshold (in years)")

# Fit generalized Pareto model with pooled data above 108
fit_gpd <- fit_elife(
  arguments = ew,
  family = "gp",
  thresh = tyr,
  export = TRUE
)

# Fit generalized Pareto model using life table approach
# Poisson model with counts and *multiplicative* offset
pois_ll_gpd <- function(par, formula, offset, data) {
  mod <- model.frame(formula, data = data)
  y <- model.response(mod)
  X <- model.matrix(formula, data = data)
  stopifnot(length(par) == ncol(X))
  pred <- as.numeric(X %*% par)
  if (any(pred < 0)) {
    return(-Inf)
  }
  sum(dpois(x = y, lambda = offset / pred, log = TRUE))
}
# Obtain coefficients estimates for the generalized Pareto model
pois_dat_gpd <- pois_dat
if (tyr > tyr_smooth) {
  pois_dat_gpd <- pois_dat |>
    dplyr::filter(time > (thresh - thresh_smooth)) |>
    dplyr::mutate(
      time = time - (thresh - thresh_smooth),
      midbin = midbin - min(midbin)
    )
}
# Optimize the log likelihood
fit_pois <- optim(
  par = c(500, -0.05),
  fn = pois_ll_gpd,
  data = pois_dat_gpd,
  offset = pois_dat_gpd$nrisk,
  formula = formula(nfail ~ time),
  method = "Nelder",
  control = list(fnscale = -1)
)
# Standardized estimates (yearly scale) are on par
pars_diff <- fit_pois$par / c(365.25, 1) - coef(fit_gpd)

# Fit exponential distribution
fit_exp <- fit_elife(
  arguments = ew,
  family = "exp",
  thresh = tyr,
  export = TRUE
)

# Compare the two nested models
anova(fit_exp, fit_gpd)


# Similar test, but with a semiparametric (smooth) alternative
smoothhazard108 <- mgcv::gam(
  nfail ~ offset(log(nrisk)) + s(midbin, bs = "ps", k = 100),
  family = poisson,
  data = pois_dat_gpd,
  method = "REML"
)
# Compare against constant model (no smooth - constant hazard)
anova(smoothhazard108)

# Generalized cross validation
smoothhazard108_GCV <- mgcv::gam(
  nfail ~ offset(log(nrisk)) + s(midbin, bs = "ps"),
  family = poisson,
  data = pois_dat_gpd
)
smooth_conf108 <- gratia:::confint.gam(
  smoothhazard108_GCV,
  parm = "s(midbin)",
  type = "simultaneous",
  trans = exp
)
# Plot of smooth - note that the dashed lines are 1 std. err. bands, transformed
plot(smoothhazard108_GCV,
  trans = exp,
  rug = FALSE,
  ylim = c(0, 5)
)
# Compare against constant model (no smooth - constant hazard)
anova(smoothhazard108_GCV)
# Probability of surviving an extra year
exp(-(exp(coef(glm(
  nfail ~ offset(log(nrisk)),
  family = poisson,
  data = pois_dat_gpd
))) * 365.25)) # mean-zero constraint on smooth params
exp(-1 / coef(fit_exp))

# Compute the generalized Pareto hazard over the grid for the plot
haz_age_seq <- seq(tyr, 115, length.out = 100)
haz_gpd <- helife(
  x = haz_age_seq - tyr,
  scale = coef(fit_gpd)["scale"],
  shape = coef(fit_gpd)["shape"],
  family = "gp"
)


####################################################
###       Lexis diagram and sampling frame       ###
####################################################

if (generateFigure) {
  # Plot a Lexis diagram manually using ggplot2
  startDate <- lubridate::ymd("1995-01-01")
  start105 <- lubridate::ymd("1999-12-28")
  end105 <- lubridate::ymd("2015-01-01")
  df <- englandwales |>
    filter(ndays >= thresh, ddate > startDate) |>
    mutate(
      entry = pmax(startDate, bdate + pmax(ltrunc1, ceiling(365.25 * tyr))),
      exit = pmin(bdate + rtrunc1, ddate),
      ageend = pmin(rtrunc1, ndays) / 365.25,
      agestart = pmax((startDate - bdate) / 365.25, pmax(ltrunc1 / 365.25, tyr))
    )
  df2 <- englandwales |>
    filter(
      ndays >= round(110 * 365.25),
      !is.na(ltrunc2)
    ) |>
    mutate(
      entry = bdate + round(110 * 365.25),
      exit = bdate + ndays,
      ageend = ndays / 365.25,
      agestart = 110
    )
  g3 <- ggplot(df) +
    geom_segment(
      aes(
        x = entry,
        xend = exit,
        y = agestart,
        yend = ageend
      ),
      alpha = 0.3
    ) +
    geom_segment(
      data = df2,
      mapping = aes(
        x = entry,
        xend = exit,
        y = agestart,
        yend = ageend
      ),
      alpha = 0.3
    ) +
    geom_point(
      data = df |> filter(ageend == (ndays / 365.25)),
      mapping = aes(x = exit, y = ageend),
      size = 0.5, alpha = 0.3
    ) +
    geom_point(
      data = df2,
      mapping = aes(x = exit, y = ageend),
      size = 0.5, alpha = 0.3
    ) +
    # geom_rug(aes(x = entry), sides = "t", alpha = 0.5, length = grid::unit(2, "mm")) +
    scale_x_date(
      date_breaks = "3 year", minor_breaks = "1 year",
      labels = date_format("%Y"),
      limits = c(
        startDate,
        lubridate::ymd("2021-01-01")
      ),
      breaks = seq.Date(
        from = startDate,
        to = lubridate::ymd("2021-01-01"), by = "2 year"
      ),
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    scale_y_continuous(
      limits = c(tyr, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = "calendar time", y = "age (in years)")
  # Plot showing the effects of truncation - only for common time window
  # Ordered by upper bound, since the sampling
  # window is fixed but the lower bound is the maximum of the
  # threshold and the lower truncation limit
  g4 <- englandwales |>
    filter(
      ndays > thresh,
      ddate >= start105,
      ddate <= end105
    ) |>
    mutate(
      ltrunc = pmax(start105 - bdate, thresh) / 365.25,
      rtrunc = (end105 - bdate) / 365.25,
      time = ndays / 365.25,
      entry = bdate + thresh,
      .keep = "none"
    ) |>
    arrange(rtrunc) |>
    mutate(id = seq_len(dplyr::n())) |>
    ggplot() +
    geom_point(mapping = aes(x = entry, y = time), shape = 4, alpha = 0.5) +
    geom_line(mapping = aes(x = entry, y = ltrunc)) +
    geom_line(mapping = aes(x = entry, y = rtrunc)) +
    scale_x_date(
      date_breaks = "2 year", minor_breaks = "1 year",
      labels = date_format("%Y"),
      limits = c(lubridate::ymd("1994-01-01"), end105),
      breaks = seq.Date(
        from = lubridate::ymd("1994-01-01"),
        to = end105,
        by = "3 years"
      ),
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    scale_y_continuous(
      limits = c(tyr, 114), minor_breaks = seq(tyr, 114, by = 0.5),
      breaks = seq(tyr, 114, by = 1L),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = paste("date of", tyr, "anniversary"),
         y = "age at death (in years)")

  g3 + g4
  ggsave("../LaTeX/figures/EnglandWales_fig1.pdf",
    width = 8,
    height = 4
  )
}


####################################################
###          Goodness of fit testing             ###
####################################################
# We bin the observations into age bands accounting
# for the truncation. A simple solution is to use
# simulate multiple age for each person and compute
# the average number of death per category.
B <- 1000L
samp_age <- matrix(nrow = fit_gpd$nexc, ncol = B)
set.seed(202401)
for (i in seq_along(fit_gpd$time)) {
  # This analysis is conditional on bounds
  samp_age[i, ] <- longevity::samp_elife(
    family = "gp",
    type2 = "ditrunc",
    n = B,
    scale = coef(fit_gpd)["scale"],
    shape = coef(fit_gpd)["shape"],
    lower = with(fit_gpd, ltrunc)[i, , drop = FALSE],
    upper = with(fit_gpd, rtrunc)[i, , drop = FALSE]
  )
}
# Observed counts (in half year increments)
rfac <- 2 # rounding factor (e.g., 2 for half year)
rounded_obs <- floor(rfac * fit_gpd$time) / rfac
breaks <- c(seq(0, 113 - fit_gpd$thresh, by = 1 / rfac), Inf)
obs_bins <- table(cut(rounded_obs, breaks = breaks, right = FALSE))
rounded_exp <- floor(rfac * samp_age) / rfac
exp_bins <- table(cut(rounded_exp, breaks = breaks, right = FALSE)) / B
K <- length(exp_bins)
isTRUE(all(exp_bins > 5))
# Compute test statistic and p-value
chisq_stat <- sum((obs_bins - exp_bins)^2 / exp_bins)
pchisq(chisq_stat, df = K - 3L, lower.tail = FALSE)

# When do people die after their birthday? There is
# seemingly a pattern: people die just after they
# were celebrated, and less as times goes...
# Might be due to confounding with time of year
# (both for birth months)
yday <- pmax(0, as.integer(
  with(englandwales, ddate - bdate - floor(365.25 * ageyear))
))
ggplot(
  data = data.frame(yday),
  mapping = aes(x = yday)
) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
    bins = 12,
    center = 16.5
  ) +
  labs(
    x = "days after birthday",
    y = "proportion"
  )


# We can recycle the life table approach for the hazard
# We use larger bins this time to avoid inaccurate large
# sample approximation (expected counts larger than 0.5)



# Bin data to account for truncation
bins <- seq(from = thresh, to = 113 * 365.25, by = 0.5 * 365.25)
# Number at risk (conditional on NOT having failed)
nrisk <- sapply(bins, function(bin) {
  with(
    englandwales,
    sum(ltrunc1 <= bin & bin < rtrunc1 & ndays >= bin) +
      sum(ltrunc2 <= bin & bin < rtrunc2 & ndays >= bin, na.rm = TRUE)
  )
})
# Number of failures
bin_death <- (englandwales$ndays[englandwales$ndays > thresh] - thresh)
# Scale time to start at threshold
nfail <- as.integer(table(
  cut(
    x = bin_death,
    breaks = c(bins, Inf) - thresh
  )
))

pbins <- mev::pgp(
  q = c(bins, Inf) / 365.25,
  loc = tyr,
  scale = coef(fit_gpd)["scale"],
  shape = coef(fit_gpd)["shape"]
)
exp_counts <- nrisk * (diff(pbins) / (1 - pbins[-length(pbins)]))
isTRUE(all(exp_counts > 5))
min(exp_counts)
# Pearson's chi-square statistic
Pchisq <- sum((exp_counts - nfail)^2 / exp_counts)
# P-value
pchisq(Pchisq, df = length(exp_counts) - 3L, lower.tail = FALSE)
# Note that this is (overly?) sensitivity to the cutoffs and the
# number of bins, but overall little evidence against the model



####################################################
###       Bayesian model and extrapolation       ###
####################################################

# Define the log posterior for the generalized Pareto model
gp_lpost <- function(pars, thresh = tyr) {
  -nll_elife(
    par = pars,
    time = ew$time,
    ltrunc = ew$ltrunc,
    rtrunc = ew$rtrunc,
    family = "gp",
    thresh = thresh
  ) +
    revdbayes::gp_mdi(pars) # MDI prior, truncated above xi > -1
}
# Use ratio-of-uniform algorithm to sample independent
# draws from the posterior
nsim <- 1e3L
# We repeat this for different thresholds
post_samp <- list()
th_seq <- seq(105, 110, by = 0.5)
for (i in seq_along(th_seq)) {
  set.seed(202401)
  # Estimate mode manually
  mode <- optim(
    par = c(1.3, -0.05),
    fn = gp_lpost,
    thresh = th_seq[i],
    control = list(fnscale = -1),
    hessian = TRUE
  )
  # Draw samples from the posterior
  post_samp[[i]] <- rust::ru(
    logf = gp_lpost,
    thresh = th_seq[i],
    n = nsim,
    d = 2,
    mode = mode$par,
    lower = c(0, -0.5),
    rotate = TRUE,
    var_names = c("scale", "shape")
  )$sim_vals
}
save(post_samp, file = "01-EnglandWales-postsamp.RData")

# Extract credible intervals for the shape parameter
cred_interv <- t(sapply(post_samp, function(x) {
  quantile(x[, 2], probs = c(0.025, 0.5, 0.975))
}))

# Compute posterior of endpoint parameter
post_endpoint <- sapply(seq_along(th_seq), function(i) {
  th_seq[i] + ifelse(post_samp[[i]][, 2] >= 0,
    Inf,
    -post_samp[[1]][, 1] / post_samp[[1]][, 2]
  )
})
# Posterior probability that endpoint exceeds 130 years
rbind(
  thresh = th_seq,
  median_endpoint = apply(post_endpoint, 2, function(x) {
    mean(x >= 130)
  })
)

# Posterior of 99% of exceedances above 110
# This calculation leverages threshold stability
post99 <- sapply(seq_along(th_seq[th_seq <= 110]), function(i) {
  revdbayes::qgp(0.99,
    loc = 110,
    scale = (110 - th_seq[i]) * post_samp[[i]][, 2] + post_samp[[i]][, 1],
    shape = post_samp[[i]][, 2]
  )
})


## Quantile-quantile plot
# 'longevity' can do this, but does not include uncertainty measures

if (generateFigure) {
  # Create plot of the hazard
  g1 <- ggplot(data = ci105_GCV) +
    geom_line(mapping = aes(
      x = time,
      y = est
    )) +
    geom_line(mapping = aes(x = time, y = lower), linetype = "dashed") +
    geom_line(mapping = aes(x = time, y = upper), linetype = "dashed") +
    geom_line(
      data = data.frame(
        x = haz_age_seq,
        y = haz_gpd
      ),
      mapping = aes(x = x, y = y),
      linewidth = 1.2,
      col = "grey"
    ) +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(
      limits = c(tyr_smooth, 115),
      breaks = seq(tyr_smooth, 115, by = 2L)
    ) +
    labs(
      x = "age (in years)",
      y = "hazard (in years)"
    ) +
    theme_classic()

  # This position is used for PP plots, and whichever plot which maps back to common scale
  source("03-helper_functions.R")
  g2 <- qqplot(
    object = fit_gpd,
    pars = post_samp[[which(th_seq == fit_gpd$thresh)]]
  )

  g1 + g2
  ggsave(
    "../LaTeX/figures/EnglandWales_fig2.pdf",
    width = 8,
    height = 4
  )

  # Parameter stability plot
  # Compare profile-likelihood and credible intervals
  g7 <- ggplot() +
    geom_linerange(
      data = data.frame(
        thresh = th_seq,
        ymin = cred_interv[, 1],
        y = cred_interv[, 2],
        ymax = cred_interv[, 3]
      ),
      mapping = aes(x = thresh, y = y, ymin = ymin, ymax = ymax)
    ) +
    # Get frequentist estimates (profile), jittered
    geom_linerange(
      data = as.data.frame(cbind(thresh = tstab$thresh + 0.05, tstab$shape)),
      mapping = aes(
        x = thresh + 0.05,
        y = estimate,
        ymin = lower,
        ymax = upper
      ),
      col = "grey"
    ) +
    labs(y = "shape", x = "threshold (in years)")

  g8 <- ggplot(
    data = data.frame(
      endpoint = c(post99),
      thresh = factor(rep(th_seq[1:ncol(post99)], each = nsim))
    ),
    mapping = aes(x = thresh, y = endpoint)
  ) +
    ggdist::stat_slabinterval() +
    scale_y_continuous(limits = c(114, 120)) +
    scale_x_discrete(labels = c(rbind(105:110, rep("", 6)))[-12]) +
    labs(
      x = "threshold (in years)",
      # subtitle = "0.99 quantile of exceedances above 110 years",
      y = "age (in years)"
    )
  g7 + g8
  ggsave("../LaTeX/figures/EnglandWales_fig3.pdf",
    width = 8,
    height = 4
  )
}


####################################################
###     Differences between cohorts/sex          ###
####################################################

# Recent Demographic Research paper fit a Gompertz
# distribution with a scale that varies according to
# sex and with a linear time for birth year

test_gp <- longevity::test_elife(
  time = ew$time,
  ltrunc = ew$ltrunc,
  rtrunc = ew$rtrunc,
  family = "gp",
  thresh = tyr,
  covariate = ew$gender
)
test_gomp <- longevity::test_elife(
  time = ew$time,
  ltrunc = ew$ltrunc,
  rtrunc = ew$rtrunc,
  family = "gomp",
  thresh = tyr,
  covariate = ew$gender
)

## Smoothing of mortality
# We fit a model to supercentenarians and use the
ew_bounds <- longevity::idlmetadata |> filter(country == "EAW", group == "110+")
thresh110 <- floor(110 * 365.25)
supercent <- englandwales |>
  dplyr::filter(ndays > thresh110) |>
  dplyr::mutate(
    ltrunc = as.numeric(pmax(thresh110, ew_bounds$ldate - bdate) - thresh110),
    rtrunc = as.numeric(ew_bounds$rdate - bdate - thresh110),
    age = ndays - thresh110
  )
sDate <- seq.Date(
  from = ymd("1990-01-01"),
  to = ymd("2020-12-31"),
  by = "month"
)


# Pick bandwidth for smoothing
bandwidth <- 365.25 * 5 # number of days

npmle <- matrix(ncol = 5L, nrow = length(sDate))
par_exp <- matrix(ncol = 3L, nrow = length(sDate))
conv <- rep(NA, length.out = length(sDate))
colnames(npmle) <- c("nobs", "Q1", "median", "Q3", "mean")
for (i in seq_along(sDate)) {
  weights <- gaussian_weight(
    x = supercent$ddate,
    x0 = sDate[i],
    bandwidth = bandwidth
  )
  # Truncate small weights to speed up fitting
  weights <- ifelse(weights < 0.01, 0, weights)
  # Compute effective number of observations
  nobs <- sum(weights)
  sdat <- supercent[weights > 0, ]
  npmle[i, "nobs"] <- nobs
  if (nobs > 10) {
    np <- suppressWarnings(np_elife(
      time = sdat$age / 365.25,
      thresh = 0,
      ltrunc = sdat$ltrunc / 365.25,
      rtrunc = sdat$rtrunc / 365.25,
      event = rep(1, nrow(sdat)),
      weights = weights[weights > 0], # Remove 0 weights to speed up routine
      method = "em", # Use CPP implementation
      maxiter = 1e4 # Reduce max number of iteration
    ))
    # Convergence errors can be spurious with large databases
    conv[i] <- np$abstol
    # Compute quartiles and restricted mean
    npmle[i, -1] <- summary(np)[c(2, 4, 5, 3)]
    # Fit an exponential model
    exp_fit <- fit_elife(
      time = sdat$age / 365.25,
      thresh = 0,
      ltrunc = sdat$ltrunc / 365.25,
      rtrunc = sdat$rtrunc / 365.25,
      event = rep(1, nrow(sdat)),
      weights = weights[weights > 0],
      family = "exp"
    )
    prof_exp <- longevity::prof_exp_scale(
      mle = exp_fit,
      time = sdat$age / 365.25,
      thresh = 0,
      ltrunc = sdat$ltrunc / 365.25,
      rtrunc = sdat$rtrunc / 365.25,
      event = rep(1, nrow(sdat)),
      weights = weights[weights > 0]
    )
    par_exp[i, ] <- prof_exp
  }
}

# Plot quartiles over time along with mean,
# estimated using the exponential model
g9 <- ggplot(data = data.frame(
  date = sDate,
  y = 110 + par_exp[, 1],
  lower = 110 + par_exp[, 2],
  upper = 110 + par_exp[, 3]
)) +
  geom_line(aes(x = date, y = y), linewidth = 1.5) +
  geom_line(aes(x = date, y = lower), linetype = "dashed") +
  geom_line(aes(x = date, y = upper), linetype = "dashed") +
  geom_line(
    data = data.frame(
      x = sDate,
      y = 110 + npmle[, "median"]
    ),
    aes(x = x, y = y),
    col = "grey"
  ) +
  geom_line(
    data = data.frame(
      x = sDate,
      y = 110 + npmle[, "Q1"]
    ),
    aes(x = x, y = y),
    col = "grey"
  ) +
  geom_line(
    data = data.frame(
      x = sDate,
      y = 110 + npmle[, "Q3"]
    ),
    aes(x = x, y = y),
    col = "grey"
  ) +
  geom_point(
    data = supercent |> dplyr::filter(ddate > lubridate::ymd("1990-01-01")),
    mapping = aes(x = ddate, y = ndays / 365.25),
    alpha = 0.2,
    shape = 4
  ) +
  scale_y_continuous(
    limits = c(110, 114), expand = c(0, 0),
    oob = scales::oob_keep
  ) +
  scale_x_date(expand = c(0, 0)) +
  labs(
    y = "age (in year)",
    x = "death date"
  ) +
  theme_classic()
g9
ggsave(
  "../LaTeX/figures/EnglandWales_fig4.pdf",
  width = 4,
  height = 4
)
