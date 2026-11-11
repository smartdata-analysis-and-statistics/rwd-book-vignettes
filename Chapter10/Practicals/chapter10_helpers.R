# Helpers for the comparative-effectiveness practicals.
# The functions deliberately generate teaching data; they are not a substitute
# for a prespecified analysis of real clinical data.

generate_ce_data <- function(n = 3000, seed = 1234, heterogeneous = FALSE,
                             quadratic_cost_age = FALSE) {
  set.seed(seed)

  # Transform normal draws to a truncated normal instead of clipping them.
  # Clipping creates artificial piles at exactly ages 18 and 64. Reusing the
  # normal draws also keeps the random-number stream for other variables fixed.
  age_draw <- rnorm(n, 46, 10)
  age_probability <- pnorm(age_draw, mean = 46, sd = 10)
  age <- qnorm(
    pnorm(18, mean = 46, sd = 10) +
      age_probability * (pnorm(64, mean = 46, sd = 10) -
                           pnorm(18, mean = 46, sd = 10)),
    mean = 46, sd = 10
  )
  gender <- factor(rbinom(n, 1, 0.75), labels = c("Male", "Female"))
  prevDMTefficacy <- factor(
    sample(c("None", "Low_efficacy", "Medium_high_efficacy"), n,
           replace = TRUE, prob = c(0.35, 0.08, 0.34)),
    levels = c("None", "Low_efficacy", "Medium_high_efficacy")
  )
  if (quadratic_cost_age) {
    # Optional teaching variant: add a U-shaped age association while keeping
    # generated costs within the range of the original seeded cost draw.
    legacy_cost_range <- c(5.319158, 13.16082)
    cost_mean <- 9 + 0.009 * ((age - 41)^2 - 130)
    logPremedicalcost <- pmin(
      pmax(rnorm(n, cost_mean, 1.15), legacy_cost_range[1]),
      legacy_cost_range[2]
    )
  } else {
    logPremedicalcost <- rnorm(n, 9, 1.15)
  }
  prerelapseNum <- rpois(n, 0.43)
  numSymptoms <- pmin(rpois(n, 1.2), 4)
  years <- runif(n, 0.5, 3)

  # Treatment is deliberately confounded with baseline characteristics.
  p_dmf <- plogis(0.8 - 0.06 * (age - 45) +
    0.35 * (gender == "Female") - 0.25 * (logPremedicalcost - 9))
  treatment <- factor(rbinom(n, 1, p_dmf), labels = c("TERI", "DMF"))

  iscore <- as.numeric(scale(age)) - 0.35 * as.numeric(gender == "Female") +
    0.25 * logPremedicalcost + rnorm(n, 0, 0.6)
  treatment_log_effect <- if (heterogeneous) {
    0.18 + 0.28 * as.numeric(scale(iscore))
  } else {
    rep(0.25, n)
  }

  log_lambda <- 0.45 + 0.015 * age + 0.20 * (gender == "Male") +
    0.12 * logPremedicalcost + 0.18 * prerelapseNum +
    0.08 * numSymptoms + treatment_log_effect * (treatment == "DMF") +
    log(years)
  y <- rpois(n, exp(log_lambda))

  data.frame(
    y, treatment, age, gender, prevDMTefficacy, logPremedicalcost,
    prerelapseNum, numSymptoms, years, Iscore = iscore
  )
}

get_missdata <- function(data,
                         scenario = c("MAR", "MAR_COVARIATE", "MART", "MARTY", "MNAR"),
                         seed = 1234) {
  scenario <- match.arg(scenario)
  set.seed(seed)
  dat <- data
  risk <- switch(
    scenario,
    # Keep the original gender-only MAR pattern for existing examples.
    MAR = plogis(-2.3 + 0.2 * (dat$gender == "Female")),
    # The covariate MAR pattern depends on fully observed predictors used by
    # the imputation model. It excludes the incomplete variables themselves.
    MAR_COVARIATE = plogis(-2.3 +
      0.2 * (dat$gender == "Female") +
      0.35 * (dat$treatment == "TERI") +
      0.25 * as.numeric(scale(dat$y)) +
      0.2 * as.numeric(scale(dat$numSymptoms)) +
      0.15 * as.numeric(scale(dat$years))),
    MART = plogis(-2.3 + 0.03 * (dat$age - 45) + 0.2 * (dat$gender == "Female") +
      0.45 * (dat$treatment == "TERI")),
    MARTY = plogis(-2.5 + 0.03 * (dat$age - 45) + 0.2 * (dat$gender == "Female") +
      0.45 * (dat$treatment == "TERI") + 0.02 * dat$y),
    MNAR = plogis(-2.3 + 0.03 * (dat$age - 45) + 0.25 * dat$logPremedicalcost)
  )
  dat$logPremedicalcost[runif(nrow(dat)) < risk] <- NA
  dat$age[runif(nrow(dat)) < risk * 0.25] <- NA
  dat$prerelapseNum[runif(nrow(dat)) < risk * 0.25] <- NA
  dat$Iscore <- NULL
  dat
}

# Generate missingness in one target covariate while retaining the complete
# data in `data`. This version is designed for teaching: learners can change
# the scenario and missingness rate, then compare each analysis with the known
# complete-data benchmark.
get_missdata_one <- function(data,
                             scenario = c("MCAR", "MAR", "MART", "MARTY", "MNAR"),
                             target = "logPremedicalcost",
                             missing_rate = 0.40,
                             seed = 1234) {
  scenario <- match.arg(scenario)
  stopifnot(target %in% names(data), missing_rate > 0, missing_rate < 1)

  set.seed(seed)
  dat <- data
  target_true <- data[[target]]
  age_z <- as.numeric(scale(data$age))
  female <- as.numeric(data$gender == "Female")
  teri <- as.numeric(data$treatment == "TERI")
  y_z <- as.numeric(scale(data$y))
  target_z <- as.numeric(scale(target_true))

  linear_predictor <- switch(
    scenario,
    MCAR = rep(0, nrow(data)),
    MAR = 0.50 * age_z + 0.35 * female,
    MART = 0.50 * age_z + 0.35 * female + 0.90 * teri,
    MARTY = 0.50 * age_z + 0.35 * female + 0.90 * teri + 0.35 * y_z,
    MNAR = 0.50 * age_z + 0.35 * female + 0.90 * target_z
  )

  # Calibrate the intercept so every scenario has approximately the requested
  # missingness rate despite having a different linear predictor.
  intercept <- uniroot(
    function(a) mean(plogis(a + linear_predictor)) - missing_rate,
    interval = c(-30, 30)
  )$root
  p_missing <- plogis(intercept + linear_predictor)
  missing <- rbinom(nrow(data), size = 1, prob = p_missing)
  dat[[target]][missing == 1] <- NA
  attr(dat, "missingness_truth") <- list(
    scenario = scenario,
    target = target,
    missing_rate = mean(missing),
    observed_indicator = 1 - missing
  )
  dat
}

# Teaching extension: generate missingness in both a confounder block and the
# outcome. The complete data are used only to create a known simulation truth;
# the returned data contain the masked values and the two observation flags.
# R_X = 1 means that all variables in x_vars are observed; R_Y = 1 means that
# the outcome is observed. Mechanisms are defined for each missingness process,
# rather than assigning one label to the whole data frame.
get_joint_missing_data <- function(data,
                                   scenario = c("MIXED_MNAR", "MAR_BOTH",
                                                "MART_BOTH", "MNAR_OUTCOME"),
                                   x_vars = c("age", "gender", "prevDMTefficacy",
                                               "logPremedicalcost", "prerelapseNum"),
                                   outcome = "y",
                                   x_missing_rates = c(
                                     age = 0.12, gender = 0.02,
                                     prevDMTefficacy = 0.23,
                                     logPremedicalcost = 0.42,
                                     prerelapseNum = 0.13
                                   ),
                                   y_missing_rate = 0.20,
                                   seed = 1234) {
  scenario <- match.arg(scenario)
  stopifnot(all(c(x_vars, outcome, "treatment", "years") %in% names(data)))
  if (length(x_missing_rates) == 1) {
    x_missing_rates <- setNames(rep(x_missing_rates, length(x_vars)), x_vars)
  }
  stopifnot(all(x_vars %in% names(x_missing_rates)),
            all(x_missing_rates[x_vars] > 0),
            all(x_missing_rates[x_vars] < 1),
            y_missing_rate > 0, y_missing_rate < 1)

  set.seed(seed)
  # The published generator returns a data.table; use a plain data.frame here
  # so the helper behaves consistently for both data.table and data.frame input.
  data <- as.data.frame(data)
  dat <- data
  n <- nrow(data)
  z <- function(x) {
    if (is.factor(x)) x <- as.numeric(as.character(x))
    as.numeric(scale(x))
  }
  age_z <- z(data$age)
  female <- as.numeric(data$gender == "Female")
  iscore_z <- z(data$Iscore)
  dmf <- as.numeric(data$treatment == "DMF")
  years_z <- z(data$years)
  calibrate <- function(lp, target) {
    uniroot(function(a) mean(plogis(a + lp)) - target,
            interval = c(-30, 30))$root
  }

  # Each confounder has its own observation process. This creates partial
  # missingness patterns instead of deleting the whole confounder block at once.
  R_x <- lapply(x_vars, function(v) {
    # Iscore, treatment, and years are complete in the simulated data. The
    # branches below use only the predictors specified for this variable and
    # avoid making a variable's own missingness depend on its hidden value.
    if (scenario == "MIXED_MNAR" && v %in% c("age", "gender")) {
      # Age and gender are MCAR: no patient-level predictors.
      lp_x <- rep(0, n)
    } else if (scenario == "MIXED_MNAR" && v == "logPremedicalcost") {
      # The cost process depends on the true baseline age and gender, treatment,
      # and years, but not on true cost itself.
      lp_x <- 0.55 * age_z + 0.35 * female + 0.85 * dmf + 0.25 * years_z
    } else {
      if (scenario == "MIXED_MNAR" &&
          v %in% c("prevDMTefficacy", "prerelapseNum")) {
        # These two variables are missing only as a function of their true
        # baseline Iscore and treatment.
        lp_x <- 0.45 * iscore_z + 0.85 * dmf
      } else {
        lp_x <- 0.45 * iscore_z + 0.25 * years_z
        if (v == "logPremedicalcost") lp_x <- lp_x + 0.25 * years_z
        if (v == "prerelapseNum") lp_x <- lp_x + 0.20 * iscore_z
        if (scenario %in% c("MART_BOTH", "MIXED_MNAR")) {
          lp_x <- lp_x + 0.85 * dmf
        }
      }
    }
    p_x <- plogis(calibrate(lp_x, x_missing_rates[[v]]) + lp_x)
    rbinom(n, 1, 1 - p_x)
  })
  names(R_x) <- paste0("R_", x_vars)
  R_X <- as.integer(rowSums(as.data.frame(R_x)) == length(x_vars))

  lp_y <- 0.30 * iscore_z + 0.65 * dmf + 0.20 * years_z
  if (scenario %in% c("MNAR_OUTCOME", "MIXED_MNAR")) {
    lp_y <- lp_y + 0.45 * z(data[[outcome]])
  }
  p_y <- plogis(calibrate(lp_y, y_missing_rate) + lp_y)
  R_Y <- rbinom(n, 1, 1 - p_y)

  dat[x_vars] <- Map(function(x, r) {
    x[r == 0] <- NA
    x
  }, dat[x_vars], R_x)
  dat[names(R_x)] <- R_x
  dat[[outcome]][R_Y == 0] <- NA
  dat$R_X <- R_X
  dat$R_Y <- R_Y
  attr(dat, "missingness_truth") <- list(
    scenario = scenario,
    x_vars = x_vars,
    outcome = outcome,
    x_missing_rates = vapply(R_x, function(r) mean(r == 0), numeric(1)),
    x_block_missing_rate = mean(R_X == 0),
    y_missing_rate = mean(R_Y == 0),
    R_X = R_X,
    R_Y = R_Y
  )
  dat
}
