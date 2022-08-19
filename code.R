library(INLA)
library(sp)


#' Calculate likelihood weights.
calc_likelihood_weights <- function(years, weight.rate) {
  stopifnot(years == sort(years))
  w <- weight.rate * (1 - weight.rate) ^ (max(years) - years)
  w <- w / max(w)
  if (max(w) / min(w) > 1e5)
    warning("Weight rate too large, may cause instability.")
  w
}

#' Calculate row summaries. Data frame containing row means, medians, standard deviations, 2.5th and 97.5th percentiles
calculate_row_summaries <- function(m) {
  qs <- apply(m, 1, quantile, probs = c(0.025, 0.5, 0.975))
  data.frame(
    mean = rowMeans(m), 
    median = qs[2, ],
    sd = apply(m, 1, sd),
    lb = qs[1, ],
    ub = qs[3, ]
  )
}

#' A convenience function to extract the mortality matrix from a model fit.
fitted_rates_matrix <- function(model.fit) UseMethod("fitted_rates_matrix")

#' @export
fitted_rates_matrix.inla <- function(model.fit) inla_fitted_matrix(model.fit, "mean")

#' @export
fitted_rates_matrix.lc <- function(model.fit) model.fit$rates

#' A convenience function to extract the mortality matrix from a model fit.
#' @param model.fit A model fitted by INLA.
#' @param variable Variable to be extracted from the fitted values matrix; mean, median, etc.
inla_fitted_matrix <- function(model.fit, variable) {
  v <- model.fit$summary.fitted.values[[variable]]
  v <- v[order(model.fit$.args$data$year, model.fit$.args$data$age)]
  m <- matrix(v, nrow = length(unique(model.fit$.args$data$age)),
              ncol = length(unique(model.fit$.args$data$year)))
  colnames(m) <- unique(model.fit$.args$data$year)
  rownames(m) <- unique(model.fit$.args$data$age)
  m
}

# Kanisto Thatcher method
kt_extension <- function(lx70) {
  mux <- (log(lx70[-nrow(lx70), , drop = FALSE]) - log(lx70)[-1, , drop = FALSE]) / 5
  
  lx70 <- rbind(
    lx70[rep(seq(3), each = 5), , drop = FALSE] *
      exp(-seq(0, 4) * mux[rep(seq(3), each = 5), , drop = FALSE]),
    lx70[4, , drop = FALSE]
  )
  
  dx70 <- rbind(lx70[-nrow(lx70), , drop = FALSE] - lx70[-1, , drop = FALSE],
                lx70[nrow(lx70), , drop = FALSE]
  )
  qx70 <- dx70 / lx70
  
  logitqx70 <- log(qx70 / (1 - qx70))
  logitqx70[nrow(logitqx70), ] <- NA # Not defined for 85+
  y <- as.vector(logitqx70)
  x <- rep(70:85 + .5, length.out = length(y))
  num.yr <- length(y) / 16
  yr <- as.factor(rep(seq_len(num.yr), each = 16))
  w <- as.vector(dx70)
  if (nlevels(yr) == 1) {
    mod <- lm(y ~ x, weights = x)
    logA <- mod$coefficients[1] # intercept
    B <- mod$coefficients[2] # slope
  } else {
    mod <- lm(y ~ 0 + yr + yr:x, weights = w)
    logA <- mod$coefficients[paste0("yr", seq(num.yr))] # intercepts
    B <- mod$coefficients[paste0("yr", seq(num.yr), ":x")] # slopes
  }
  
  logitqx85 <- t(logA + outer(B, (85:129 + .5))) # predict logit qx for age >= 85
  qx85 <- exp(logitqx85) / (1 + exp(logitqx85)) # invert logit transform
  
  lx85 <- matrix(nrow = nrow(logitqx85), ncol = ncol(logitqx85))
  lx85[1, ] <- lx70[nrow(lx70), ] # last entry of vector holding lx70-85
  for (k in seq(2, nrow(lx85)))
    lx85[k, ] <- lx85[k - 1, ] * (1 - qx85[k - 1, ])
  
  dx85 <- lx85 * qx85
  
  lx70 <- rbind(lx70[-nrow(lx70), , drop = FALSE], lx85)
  dx70 <- rbind(dx70[-nrow(dx70), , drop = FALSE], dx85)
  
  yl <- 70:129 - c(rep(c(70, 75, 80), each = 5), rep(85, length(85:129))) + 0.5
  
  x5y <- seq(70, 85, 5)[findInterval(70:129, seq(70, 85, 5))]
  
  ax70 <- as.vector(t(sapply(split(seq(nrow(dx70)), x5y),
                             function(v) colSums(dx70[v, , drop = FALSE] * yl[v]) /
                               colSums(dx70[v, , drop = FALSE]))
  ))
  ax70
}

# Model averaged projections of life expectancy

maple2 <- function(deaths, population, forecast.horizon, holdout, models = maple_models(),
                  num.draws = 1000, ax = NULL, num.threads = inla.getOption("num.threads"),
                  verbose = TRUE) {
  
  holdout.cols <- ncol(deaths) + seq(-holdout + 1, 0)
  
  if (verbose) message("Fitting models for model weight calculation...")
  weight.run.fits <- maple_fit_ensemble(
    deaths = deaths[, -holdout.cols],
    population = population[, -holdout.cols],
    forecast.horizon = holdout,
    models = models,
    num.draws = num.draws,
    ax = ax[, -holdout.cols],
    num.threads = num.threads,
    verbose = verbose)
  if (verbose) message("Forecasting with individual  models...")
  forecast.run.fits <- maple_fit_ensemble(
    deaths = deaths,
    population = population,
    forecast.horizon = forecast.horizon,
    models = models,
    num.draws = num.draws,
    ax = ax,
    num.threads = num.threads,
    verbose = verbose)
  
  if (verbose) message("Calculating projection errors...")
  projection.errors <- sapply(weight.run.fits$model.fits, function(m) {
    maple_projection_error(
      deaths = deaths[, holdout.cols],
      population = population[, holdout.cols],
      ax = ax[, holdout.cols],
      fitted.rates = fitted_rates_matrix(m)[, holdout.cols]
    )
  })
  if (verbose) message("Computing model average...")
  model.weights <- maple_model_weights(projection.errors)
  
  model.draws <- floor(model.weights * floor(min(num.draws / model.weights)))
  
  draw.idx <- lapply(seq_along(model.draws), function(i) {
    sample(num.draws, model.draws[i], replace = FALSE)
  })
  
  bma.samples <- unlist(
    lapply(seq_along(forecast.run.fits$samples),
           function(i) forecast.run.fits$samples[[i]][draw.idx[[i]]]),
    recursive = FALSE)
  
  bma.sample.summaries <- maple_sample_summaries(bma.samples)
  
  bma <- structure(list(model.weights = model.weights,
                        sample.summaries = bma.sample.summaries,
                        samples = bma.samples,
                        individual.model.sample.summaries = forecast.run.fits$sample.summaries),
                   class = "maple.bma")
  if (verbose) message("Done.")
  bma
}

#' Fit an ensemble of models.

maple_fit_ensemble <- function(deaths, population, forecast.horizon, models = maple_models(),
                               num.draws = 1000, ax = NULL, num.threads = inla.getOption("num.threads"),
                               verbose = TRUE) {
  
  
  
  model.fits <- structure(lapply(seq_along(models), function(m) {
    if (verbose) message("Fitting model ", models[[m]]$name, "...")
    tryCatch(maple_fit_model(model = models[[m]],
                             deaths = deaths,
                             population = population,
                             forecast.horizon = forecast.horizon,
                             num.threads = num.threads),
             error = function(e) {
               message("Error:", e$message)
               return(e)
             }
    )
  }), 
  names = names(models)
  )
  
  samples <- summaries <- NULL
  if (num.draws > 0) {
    samples <- structure(lapply(seq_along(models), function(m) {
      if (verbose) message("Sampling draws for model ", models[[m]]$name, "...")
      tryCatch(
        maple_sample(model.fit = model.fits[[m]],
                     num.draws = num.draws,
                     ax = ax),
        error = function(e) {
          message("Error:", e$message)
          return(e)
        }
      )
    }),
    names = names(models)
    )
    
    sample.summaries <- lapply(seq_along(models), function(m) {
      if (verbose) message("Calculating summary statistics for model ",
                           models[[m]]$name, "...")
      tryCatch(
        data.frame(model = names(models)[m], 
                   maple_sample_summaries(samples[[m]])),
        error = function(e) {
          message("Error:", e$message)
          return(e)
        }
      )
    })
    sample.summaries <- do.call(rbind, sample.summaries)
  }
  list(
    sample.summaries = sample.summaries,
    samples = samples,
    model.fits = model.fits
  )
}

#' Fit a single forecasting model.

maple_fit_model <- function(model, deaths, population, forecast.horizon, num.threads) {
  UseMethod("maple_fit_model")
}

maple_fit_model.inla.model <- function(model, deaths, population, forecast.horizon,
                                       num.threads = inla.getOption("num.threads"), ...) {
  inla.dat <- prep_inla_dat(deaths = deaths,
                            population = population,
                            model = model,
                            forecast.horizon = forecast.horizon)
  lhd.wght <- NULL
  if (!is.null(model$likelihood.weight.rate)) {
    inla.setOption(enable.inla.argument.weights = TRUE)
    lhd.wght <- calc_likelihood_weights(inla.dat$year, model$likelihood.weight.rate)
  }
  fit <- inla(model$fml,
              family = "xpoisson",
              data = inla.dat,
              E = population,
              control.predictor = list(link = 1),
              control.compute = list(dic = TRUE, config = TRUE),
              weights = lhd.wght,
              num.threads = num.threads
  )
  if (!is.null(model$likelihood.weight.rate)) {
    inla.setOption(enable.inla.argument.weights = FALSE)
  }
  fit
}

maple_fit_model.lc.model <- function(model, deaths, population, forecast.horizon, ...) {
  ages <- as.numeric(rownames(deaths))
  years <- as.numeric(colnames(deaths))
  forecast.years <- ncol(deaths) + seq(forecast.horizon)
  
  rates <- deaths / population
  
  if (any(is.na(rates))) stop("Cannot run Lee-Carter model beacause of NAs in data.")
  
  logrates <- log(rates)
  alphas <- rowMeans(logrates)
  clogrates <- logrates - alphas
  clogrates.svd <- svd(clogrates)
  betas <- t(clogrates.svd$u[, seq_len(model$num.pcs)])
  gammas <- betas %*% clogrates
  # Normalise so that each beta sums to 1
  gammas <- gammas * rowSums(betas)
  betas <- betas / rowSums(betas)
  gammas.fit <- apply(gammas, 1, arima, order = c(0, 1, 0), xreg = seq_len(ncol(gammas)))
  # In-sample fits
  logrates.fit <- alphas + Reduce(`+`, lapply(seq_len(nrow(gammas)),
                                              function(i) outer(betas[i, ], gammas[i, ])))
  # Forecasts
  gammas.pred <- lapply(gammas.fit, predict,
                        newxreg = ncol(gammas) + seq(forecast.horizon),
                        n.ahead = forecast.horizon)
  gammas.pred.mean <- lapply(gammas.pred, function(x) as.numeric(x$pred))
  gammas.pred.se <- lapply(gammas.pred, function(x) as.numeric(x$se))
  logrates.pred <- alphas + Reduce(`+`, lapply(seq_along(gammas.pred.mean),
                                               function(i) outer(betas[i, ], gammas.pred.mean[[i]])))
  
  logfit <- cbind(logrates.fit, logrates.pred)
  rownames(logfit) <- ages
  colnames(logfit) <- c(years, years[length(years)] + seq(forecast.horizon))
  fit <- exp(logfit)
  
  l <- structure(list(
    logrates = logfit,
    rates = fit,
    alphas = alphas,
    betas = betas,
    gammas = gammas,
    gammas.fit = gammas.fit,
    gammas.pred = gammas.pred,
    pct.var = sum(clogrates.svd$d[seq_len(model$num.pcs)] ^ 2 / sum(clogrates.svd$d ^ 2)),
    num.pcs = model$num.pcs,
    forecast.horizon = forecast.horizon
  ), class = "lc")
  l
}

#' Calculate model weights.

maple_model_weights <- function(projection.errors) {
  exp(- abs(projection.errors)) / sum(exp(- abs(projection.errors)))
}

#' Generate a list of forecasting models.

maple_models <- function(fixed.prec = 0.001,
                         gamma.shape = 1,
                         gamma.rate = 1e-3) {
  
  hpr <- paste0('list(prec = list(prior = "loggamma", param = c(', 
                gamma.shape, ', ', gamma.rate, ')))')
  
  models <- list(
    list(
      name = "IIDAGE",
      desc = "Age-period model with iid age intercepts and slopes.",
      fml = 'deaths ~
            f(age.id1, model = "iid", hyper = hpr) +
            f(yearc1, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc1, model = "iid", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id3, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "RW1AGE",
      desc = "Age-period model with first order random walk age intercepts and slopes.",
      fml = 'deaths ~
            f(age.id1, model = "rw1", hyper = hpr) +
            f(yearc1, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc1, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id3, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "IIDAGE_RW1COH",
      desc = "Age-period-cohort model with iid age intercepts and slopes and first order random walk cohort slopes.",
      fml = 'deaths ~
            f(age.id1, model = "iid", hyper = hpr) +
            f(yearc1, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc1, model = "iid", hyper = hpr) +
            f(cohort.id1, yearc1, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id3, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "RW1AGE_RW1COH",
      desc = "Age-period-cohort model with first order random walk age intercepts and slopes and first order random walk cohort slopes.",
      fml = 'deaths ~
            f(age.id1, model = "rw1", hyper = hpr) +
            f(yearc1, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc, model = "rw1", hyper = hpr) +
            f(cohort.id1, yearc, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id3, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL10_IIDAGE",
      desc = "Age-period model with iid age intercepts and slopes, using piecewise linear slopes with knot 10 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "iid", hyper = hpr) +
            f(yearc1.10a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.10b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.10a, model = "iid", hyper = hpr) +
            f(age.id3, yearc2.10b, model = "iid", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL10_RW1AGE",
      desc = "Age-period model with first order random walk age intercepts and slopes, using piecewise linear slopes with knot 10 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "rw1", hyper = hpr) +
            f(yearc1.10a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.10b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.10a, model = "rw1", hyper = hpr) +
            f(age.id3, yearc2.10b, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL6_IIDAGE",
      desc = "Age-period model with iid age intercepts and slopes, using piecewise linear slopes with knot 6 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "iid", hyper = hpr) +
            f(yearc1.6a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.6b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.6a, model = "iid", hyper = hpr) +
            f(age.id3, yearc2.6b, model = "iid", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL6_RW1AGE",
      desc = "Age-period model with first order random walk age intercepts and slopes, using piecewise linear slopes with knot 6 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "rw1", hyper = hpr) +
            f(yearc1.6a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.6b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.6a, model = "rw1", hyper = hpr) +
            f(age.id3, yearc2.6b, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL10_IIDAGE_RW1COH",
      desc = "Age-period-cohort model with iid age intercepts and slopes and first order random walk cohort slopes, using piecewise linear slopes with knot 10 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "iid", hyper = hpr) +
            f(yearc1.10a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.10b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.10a, model = "iid", hyper = hpr) +
            f(age.id3, yearc2.10b, model = "iid", hyper = hpr) +
            f(cohort.id1, yearc3.10a, model = "rw1", hyper = hpr) +
            f(cohort.id2, yearc3.10b, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL10_RW1AGE_RW1COH",
      desc = "Age-period-cohort model with first order random walk age intercepts and slopes and first order random walk cohort slopes, using piecewise linear slopes with knot 10 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "rw1", hyper = hpr) +
            f(yearc1.10a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.10b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.10a, model = "rw1", hyper = hpr) +
            f(age.id3, yearc2.10b, model = "rw1", hyper = hpr) +
            f(cohort.id1, yearc3.10a, model = "rw1", hyper = hpr) +
            f(cohort.id2, yearc3.10b, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL6_IIDAGE_RW1COH",
      desc = "Age-period-cohort model with iid age intercepts and slopes and first order random walk cohort slopes, using piecewise linear slopes with knot 6 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "iid", hyper = hpr) +
            f(yearc1.6a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.6b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.6a, model = "iid", hyper = hpr) +
            f(age.id3, yearc2.6b, model = "iid", hyper = hpr) +
            f(cohort.id1, yearc3.6a, model = "rw1", hyper = hpr) +
            f(cohort.id2, yearc3.6b, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    ),
    list(
      name = "PWL6_RW1AGE_RW1COH",
      desc = "Age-period-cohort model with first order random walk age intercepts and slopes and first order random walk cohort slopes, using piecewise linear slopes with knot 6 years before the end of data.",
      fml = 'deaths ~
            f(age.id1, model = "rw1", hyper = hpr) +
            f(yearc1.6a, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(yearc1.6b, model = "linear", mean.linear = 0, prec.linear = fixed.prec) +
            f(age.id2, yearc2.6a, model = "rw1", hyper = hpr) +
            f(age.id3, yearc2.6b, model = "rw1", hyper = hpr) +
            f(cohort.id1, yearc3.6a, model = "rw1", hyper = hpr) +
            f(cohort.id2, yearc3.6b, model = "rw1", hyper = hpr) +
            f(year.id1, model = "rw1", hyper = hpr, group = age.id4, control.group = list(model = "exchangeable")) +
            f(epsilon.id1, model = "iid", hyper = hpr)'
    )
  )
  
  models <- lapply(models, function(m) {
    class(m) <- "inla.model"
    m$fml <- gsub("hpr", hpr, m$fml)
    m$fml <- gsub("fixed.prec", fixed.prec, m$fml)
    m$fml <- as.formula(m$fml)
    m
  })
  wl.models <- Filter(
    function(x) x$name %in% c("IIDAGE", "RW1AGE", "IIDAGE_RW1COH", "RW1AGE_RW1COH"),
    models)
  wl.models <- lapply(wl.models, function(x) {
    x$name <- paste0("WL0.05_", x$name)
    x$likelihood.weight.rate <- 0.05
    x$desc <- paste0("Weighted likelihood ", tolower(substr(x$desc, 1, 1)), substr(x$desc, 2, nchar(x$desc)))
    x
  })
  
  lc.models <- lapply(seq(5), function(n) {
    m <- list(
      name = paste0("LC_", n, "PC"),
      desc = paste0("Lee-Carter model with ", n, " principal component",
                    ifelse(n > 1, "s", ""), "."),
      num.pcs = n
    )
    class(m) <- "lc.model"
    m
  })
  l <- c(models, wl.models, lc.models)
  setNames(l, sapply(l, `[[`, "name"))
}

#' Generate a period life table.

maple_plt <- function(death.rates, ax = NULL, check.conv = FALSE, full.table = FALSE) {
  year <- rep(as.numeric(colnames(death.rates)), each = nrow(death.rates))
  age <- rep(as.numeric(rownames(death.rates)), ncol(death.rates))
  mx <- as.vector(death.rates)
  
  if (is.null(ax)) {
    ax <- ifelse(age == 0, .5, 2.5)
  } else if (ncol(ax) < ncol(death.rates)) {
    ax <- ax[, c(seq(1, ncol(ax)), rep(ncol(ax), ncol(death.rates) - ncol(ax)))]
  }
  ax <- as.vector(ax)
  
  # Remove NA rates (this is useful when eg a given year has no mortality data;
  # the removed entries are re-inserted as NAs at the end of the function call)
  na.idx <- NULL
  if (any(is.na(mx))) {
    na.idx <- which(is.na(mx))
    stopifnot(age[na.idx] == seq(0, 85, 5))
    age0 <- age
    ax0 <- ax
    mx0 <- mx
    age <- age[-na.idx]
    mx <- mx[-na.idx]
    ax <- ax[-na.idx]
  }
  
  # Replace zero rates by a small number to avoid division by zero
  mx[mx == 0] <- 1e-10
  # Probability of dying between age x and x+4
  qx <- 5 * mx / (1 + (5 - ax) * mx)
  # If probability of dying is >1, set it to "almost" 1
  qx[qx > 1] <- 0.99999
  qx[age == 85] <- 1 # by definition
  px <- 1 - qx # probability of surviving to next age group
  lx <- rep(NA, length(px))
  lx[age == 0] <- 100000
  for (k in seq(5, 85, 5))
    lx[age == k] <- lx[age == k - 5] * px[age == k - 5]
  dx <- lx * qx
  ax[age >= 70] <- kt_extension(matrix(lx[age >= 70], nrow = 4))
  
  num.iter <- 4 # Number of iterations - see Preston et al. p.47
  iter.dat <- vector("list", num.iter + 1)
  iter.dat[[1]] <- list(ax = ax, qx = qx, lx = lx, dx = dx)
  for (i in seq(num.iter)) {
    ax.new <- ax
    for (k in seq(5, 80, 5))
      ax.new[age == k] <- (-5 / 24 * dx[age == k - 5] +
                             2.5 * dx[age == k] + 5 / 24 * dx[age == k + 5]) / dx[age == k]
    ax.new[age <= 10 | age >= 70] <- ax[age <= 10 | age >= 70]
    ax <- ax.new
    qx <- 5 * mx / (1 + (5 - ax) * mx)
    qx[qx > 1] <- 0.99999
    qx[age == 85] <- 1
    px <- 1 - qx
    lx <- rep(NA, length(px))
    lx[age == 0] <- 100000
    for (k in seq(5, 85, 5))
      lx[age == k] <- lx[age == k - 5] * px[age == k - 5]
    dx <- lx * qx
    # save result of current iteration
    iter.dat[[i + 1]] <- list(ax = ax, qx = qx, lx = lx, dx = dx)
  }
  
  if (check.conv) {
    ax.iter <- sapply(iter.dat, `[[`, "ax")
    stopifnot(ax.iter[, num.iter] - ax.iter[, num.iter - 1] < 0.01)
  }
  iter.result <- iter.dat[[num.iter + 1]]
  ax <- iter.result$ax
  qx <- iter.result$qx
  lx <- iter.result$lx
  dx <- iter.result$dx
  
  Lx <- rep(NA, length(age))
  for (k in seq(0, 80, 5))
    Lx[age == k] <- 5 * lx[age == k + 5] + ax[age == k] * dx[age == k]
  Lx[age == 85] <- lx[age == 85] / mx[age == 85]
  Tx <- rep(NA, length(age))
  Tx[age == 85] <- Lx[age == 85]
  for (k in rev(seq(0, 80, 5)))
    Tx[age == k] <- Tx[age == k + 5] + Lx[age == k]
  ex <- Tx / lx
  
  # Re-insert missing values
  if(!is.null(na.idx)) {
    mx1 <- mx0
    mx1[-na.idx] <- mx
    mx <- mx1
    ax1 <- ax0
    ax1[-na.idx] <- ax
    ax <- ax1
    age <- age0
    qx1 <- rep(NA, length(mx0))
    qx1[-na.idx] <- qx
    qx <- qx1
    ex1 <- rep(NA, length(mx0))
    ex1[-na.idx] <- ex
    ex <- ex1
    Tx1 <- rep(NA, length(mx0))
    Tx1[-na.idx] <- Tx
    Tx <- Tx1
    Lx1 <- rep(NA, length(mx0))
    Lx1[-na.idx] <- Lx
    Lx <- Lx1
    lx1 <- rep(NA, length(mx0))
    lx1[-na.idx] <- lx
    lx <- lx1
  }
  if (full.table) return(data.frame(year = year, age = age, ax = ax, mx = mx, qx = qx,
                                    ex = ex, Tx = Tx, Lx = Lx, lx = lx))
  data.frame(age = age, year = year, mx = mx, qx = qx, ex = ex)
}

#' Calculate projection error.

maple_projection_error <- function(deaths, population, ax = NULL, fitted.rates) {
  
  if (!all(dim(deaths) == dim(population) & 
           dim(deaths) == dim(fitted.rates))) 
    stop("Observed deaths, population and fitted rates matrices must have equal dimensions.")
  
  if (!all(colnames(deaths) == colnames(population) & 
           colnames(deaths) == colnames(fitted.rates))) 
    stop("Years of deaths, population and fitted rates don't match.")
  
  observed.plt <- maple_plt(death.rates = deaths / population, ax = ax)
  fitted.plt <- maple_plt(death.rates = fitted.rates, ax = ax)
  
  observed.e0 <- plt_ex(observed.plt, x = 0)
  fitted.e0 <- plt_ex(fitted.plt, x = 0)
  
  mean(observed.e0 - fitted.e0, na.rm = TRUE)
}

#' Generate posterior samples of death rates and life tables.

maple_sample <- function(model.fit, num.draws, ax) {
  UseMethod("maple_sample")
}

#' @export
maple_sample.inla <- function(model.fit, num.draws, ax = NULL) {
  log.rate.draws <- inla.posterior.sample(num.draws, model.fit)
  rate.draws <- lapply(log.rate.draws, function(draw) {
    rates <- exp(draw$latent[grep("Predictor", rownames(draw$latent))])
    rates <- rates[order(model.fit$.args$data$year, model.fit$.args$data$age)]
    rates.m <- matrix(rates, ncol = length(unique(model.fit$.args$data$year)))
    rownames(rates.m) <- unique(model.fit$.args$data$age)
    colnames(rates.m) <- unique(model.fit$.args$data$year)
    rates.m
  })
  structure(lapply(rate.draws, maple_plt, ax = ax, full.table = FALSE),
            class = "inla.samples")
}

#' @export
maple_sample.lc <- function(model.fit, num.draws, ax = NULL) {
  forecast.horizon <- model.fit$forecast.horizon
  log.rate.draws <- replicate(num.draws,
                              model.fit$alphas + Reduce(`+`, lapply(seq(model.fit$num.pcs),
                                                                    function(i) outer(model.fit$betas[i, ],
                                                                                      model.fit$gammas.pred[[i]]$pred +
                                                                                        cumsum(rnorm(forecast.horizon,
                                                                                                     sd = sqrt(model.fit$gammas.fit[[i]]$sigma2)))))),
                              simplify = FALSE)
  in.sample.fit <- model.fit$rate[, seq_len(ncol(model.fit$gammas))]
  rate.draws <- lapply(log.rate.draws, function(x) {
    m <- cbind(in.sample.fit, exp(x))
    colnames(m) <- colnames(model.fit$logrates)
    m
  })
  structure(lapply(rate.draws, maple_plt, ax = ax, full.table = FALSE), 
            class = "lc.samples")
}

#' Calculate summaries for death rates, life expectancy and probability of dying.

maple_sample_summaries <- function(samples) {
  rate.samples <- do.call(cbind, lapply(samples, `[[`, "mx"))
  ex.samples <- do.call(cbind, lapply(samples, `[[`, "ex"))
  qx.samples <- do.call(cbind, lapply(samples, `[[`, "qx"))
  
  rbind(data.frame(samples[[1]][c("year", "age")], 
                   metric = "rate", 
                   calculate_row_summaries(rate.samples)),
        data.frame(samples[[1]][c("year", "age")], 
                   metric = "ex", 
                   calculate_row_summaries(ex.samples)),
        data.frame(samples[[1]][c("year", "age")], 
                   metric = "qx", 
                   calculate_row_summaries(qx.samples)))
}

#' Extract the life expectancy values for a given age from a life table.

plt_ex <- function(plt, x = 0) {
  if (!all(c("age", "year", "ex") %in% names(plt))) {
    stop("Life table must include 'age', 'year', 'ex' columns.")
  }
  setNames(plt[plt$age == x, ]$ex, plt[plt$age == x, ]$year)
}

#' Extract the probability of dying by a given age from a life table.

plt_qx <- function(plt, x = 70) { 
  if (!all(c("age", "year", "qx") %in% names(plt))) {
    stop("Life table must include 'age', 'year', 'qx' columns.")
  }
  
  calculate_pod <- function(m) {
    # m matrix of age x year qx
    1 - Reduce(`*`, as.data.frame(1 - t(m)))
  }
  
  plt1 <- plt[plt$age < x, ]
  qx <- plt1$qx[order(plt1$year, plt1$age)]
  m <- matrix(qx, ncol = length(unique(plt1$year)))
  setNames(calculate_pod(m), plt[plt$age == x, ]$year)
}

# Format data to run INLA model.
prep_inla_dat <- function(deaths, population, model, forecast.horizon) {
  dat <- data.frame(
    year = rep(as.numeric(colnames(deaths)), each = nrow(deaths)),
    age = rep(as.numeric(rownames(deaths)), ncol(deaths)),
    deaths = as.vector(deaths),
    population = as.vector(population)
  )
  dat$death.rate <- dat$deaths / dat$population
  dat$forecast.period <- FALSE
  forecast.years <- max(dat$year) + seq(forecast.horizon)
  forecast.dat <- lapply(forecast.years, function(y) {
    d <- dat[dat$year == min(dat$year), ]
    d$deaths <- d$death.rate <- NA
    d$year <- y
    d$forecast.period <- TRUE
    d
  })
  forecast.dat <- do.call(rbind, forecast.dat)
  inla.dat <- rbind(dat, forecast.dat)
  inla.dat$cohort <- inla.dat$year - inla.dat$age
  
  # Parse model formula and create variables
  fmlstr <- as.character(model$fml[3])
  inla.dat$epsilon <- seq_len(nrow(inla.dat))
  
  idx.vars <- as.list(inla.dat[c("age", "cohort", "year", "epsilon")])
  idx.dat <- lapply(seq_along(idx.vars), function(i) {
    var <- idx.vars[[i]]
    s <- paste0(names(idx.vars)[i], ".id")
    num.fx <- max(sapply(regmatches(fmlstr, gregexpr(s, fmlstr)), length))
    if (num.fx == 0) return(NULL)
    setNames(
      data.frame(replicate(num.fx, as.integer(factor(var)), simplify = FALSE)),
      paste0(s, seq(num.fx)))
  })
  idx.dat <- do.call(cbind, Filter(Negate(is.null), idx.dat))
  inla.dat <- cbind(inla.dat, idx.dat)
  
  time.vars <- unique(unlist(regmatches(fmlstr,
                                        gregexpr("yearc[0-9]*(\\.[0-9]*[ab])?", fmlstr))))
  lin.vars <- unique(unlist(regmatches(time.vars,
                                       gregexpr("^yearc[0-9]*$", time.vars))))
  pw.vars <- unique(unlist(regmatches(time.vars,
                                      gregexpr("^yearc[0-9]*\\.[0-9]*[ab]$", time.vars))))
  knot.pos <- unique(as.numeric(gsub("[ab]", "",
                                     unlist(lapply(strsplit(pw.vars, "\\."), `[[`, 2)))))
  
  inla.dat$yearc1 <- inla.dat$year - mean(inla.dat$year)
  for (var in setdiff(lin.vars, "yearc1")) inla.dat[[var]] <- inla.dat$yearc1
  
  for (kp in knot.pos) {
    vars.a <- unlist(regmatches(pw.vars,
                                gregexpr(paste0("^yearc[0-9]*\\.", kp, "a$"), pw.vars)))
    vars.b <- unlist(regmatches(pw.vars,
                                gregexpr(paste0("^yearc[0-9]*\\.", kp, "b$"), pw.vars)))
    maxyear <- max(inla.dat[!inla.dat$forecast.period, ]$year) - kp
    maxyearc <- max(inla.dat[inla.dat$year <= maxyear, ]$yearc1)
    inla.dat[vars.a] <- ifelse(inla.dat$year <= maxyear, inla.dat$yearc1, maxyearc)
    inla.dat[vars.b] <- ifelse(inla.dat$year <= maxyear, 0, inla.dat$yearc1 - maxyearc)
  }
  
  inla.dat
}

#' @export
print.maple.bma <- function(x, ...) {
  cat("Model averaged projections of mortality and life expectancy between ",
      paste(range(x$sample.summaries$year), collapse = " and "), ".\n", sep = "")
  cat("Weights of models averaged:\n")
  print(x$model.weights)
  invisible(x)
}

#' @export
print.inla.model <- function(x, ...) {
  cat("Model", x$name, ": ")
  cat(x$desc, "\n")
  cat("INLA formula:\n")
  print(x$fml)
  if (!is.null(x$likelihood.weight.rate))
    cat("Weight rate:", x$likelihood.weight.rate, "\n")
  invisible(x)
}

#' @export
print.lc.model <- function(x, ...) {
  cat("Model", x$name, ": ")
  cat(x$desc, "\n")
  invisible(x)
}


