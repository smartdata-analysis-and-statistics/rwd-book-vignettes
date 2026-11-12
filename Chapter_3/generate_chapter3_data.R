# Generate the shared cohort for the Chapter 3 propensity-score practicals.
# `generate_data()` is the simulator from Chapter 9 resources, supplied in the
# user's request. Its implementation assigns treatment from baseline
# covariates, despite the stale comment saying that treatment is randomized.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument)) {
  dirname(normalizePath(sub("^--file=", "", script_argument[[1]])))
} else {
  normalizePath(".")
}

source(file.path(
  script_directory,
  "../Chapter10/Practicals/resources/chapter 09/functions.r"
))

raw_data <- generate_data(
  n = 10000,
  # Match the seed used by Chapter 6 of the reference book.
  seed = 1854,
  beta = rep(-0.2, 5)
)

chapter3_data <- data.frame(
  age = raw_data$age,
  female = as.integer(raw_data$gender == "Female"),
  prevDMTefficacy = raw_data$prevDMTefficacy,
  premedicalcost = exp(raw_data$logPremedicalcost),
  numSymptoms = raw_data$numSymptoms,
  prerelapse_num = raw_data$prerelapseNum,
  treatment = factor(
    ifelse(raw_data$treatment == "DMF", "DMT1", "DMT0"),
    levels = c("DMT0", "DMT1")
  ),
  y = raw_data$y,
  years = raw_data$years,
  Iscore = raw_data$Iscore
)

data_directory <- file.path(script_directory, "data")
dir.create(data_directory, showWarnings = FALSE)
output_path <- file.path(data_directory, "chapter3_ps_data.rds")
saveRDS(chapter3_data, output_path)

message("Saved ", nrow(chapter3_data), " rows to ", output_path)
