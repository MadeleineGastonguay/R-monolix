# Synthetic viral-load data simulated from the model itself (2regime_model_delay.txt), with
# known true values, on a design meant to make every parameter estimable. Unlike
# tests/testthat/fixtures/synthetic_vl_data.txt (a made-up curve, for unit tests), these data
# carry information on the immune parameters, so the Fisher information matrix can be inverted.
#
# Run in x86_64 R from the repo root (needs Monolix: Simulx for the curves, Monolix for the
# check fit). Writes, next to this script:
#   synthetic_vl_data.txt      the data (id, time, vl, censor, tlr7, vac)
#   true_parameters.csv        population values used
#   true_individual_parameters.csv
# and, to dev/step8/ (not committed), a plot of the data and the check fit's log.

devtools::load_all()   # this repo's code (create_project(), theme_report(), ...)
MONOLIX <- "/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/"
HERE <- "inst/examples/synthetic_identifiable"
MODEL <- file.path(HERE, "2regime_model_delay.txt")
set.seed(2026)

# True values ####
pop <- c(lambda = 100, betas = 7, ta = 4, NPs = 10, p = 3, m = 100)
omega <- c(lambda = 0.3, betas = 0.3, ta = 0.3, NPs = 0.3, p = 0.3, m = 0.3)  # SD of log values
beta_ta_vac <- log(2)    # vaccinated animals: ta twice as large
error_a <- 0.3           # constant error on log10 VL
lod <- 1.7               # detection limit (log10)
n_per_arm <- 10
# Dense over rebound and peak (about days 14-49), then weekly to day 84 so late (vaccinated)
# rebounders also reach their peak and set point. A first design (every 2 days to day 14,
# weekly to day 56) missed most peaks: lambda_pop and m_pop had no standard error.
times <- c(0, 7, 14, 17, 21, 24, 28, 31, 35, 38, 42, 45, 49, 56, 63, 70, 77, 84)

# Animals and their parameters ####
arms <- expand.grid(tlr7 = 0:1, vac = 0:1)
animals <- data.frame(id = seq_len(n_per_arm * nrow(arms)),
                      arms[rep(seq_len(nrow(arms)), each = n_per_arm), ], row.names = NULL)
indiv <- data.frame(id = animals$id)
for (k in names(pop)) {
  effect <- if (k == "ta") beta_ta_vac * animals$vac else 0
  indiv[[k]] <- pop[[k]] * exp(effect + rnorm(nrow(animals), 0, omega[[k]]))
}

# Model curves with Simulx ####
lixoftConnectors::initializeLixoftConnectors(software = "simulx", path = MONOLIX, force = TRUE)
lixoftConnectors::newProject(modelFile = MODEL)
indiv_file <- tempfile("indiv", fileext = ".csv"); write.csv(indiv, indiv_file, row.names = FALSE)
lixoftConnectors::defineIndividualElement(name = "true_indiv", element = indiv_file)
lixoftConnectors::defineOutputElement(name = "vl_times", element = list(data = data.frame(time = times), output = "tot"))
group <- lixoftConnectors::getGroups()[[1]]$name
lixoftConnectors::setGroupElement(group = group, elements = c("true_indiv", "vl_times"))
lixoftConnectors::setSamplingMethod("keepOrder")
lixoftConnectors::setGroupSize(group, nrow(indiv))
lixoftConnectors::runSimulation()
sim <- lixoftConnectors::getSimulationResults()$res$tot
if ("original_id" %in% names(sim)) sim$id <- sim$original_id
sim <- sim[, c("id", "time", "tot")]

# Observations: error, then censoring below the detection limit ####
d <- merge(sim, animals, by = "id")
d <- d[order(d$id, d$time), ]
obs <- d$tot + rnorm(nrow(d), 0, error_a)
d$censor <- as.integer(obs <= lod)
d$vl <- round(ifelse(d$censor == 1, lod, obs), 3)
d <- d[, c("id", "time", "vl", "censor", "tlr7", "vac")]
write.table(d, file.path(HERE, "synthetic_vl_data.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(data.frame(parameter = c(paste0(names(pop), "_pop"), paste0("omega_", names(omega)), "beta_ta_vac_1", "a"),
                     value = c(pop, omega, beta_ta_vac, error_a)),
          file.path(HERE, "true_parameters.csv"), row.names = FALSE)
write.csv(merge(animals, indiv, by = "id"), file.path(HERE, "true_individual_parameters.csv"), row.names = FALSE)

cat("Animals:", nrow(animals), "; observations:", nrow(d), "; censored:", sum(d$censor),
    sprintf("(%.0f%%)", 100 * mean(d$censor)), "\n")
print(aggregate(cbind(censored = censor) ~ time, d, mean))   # share censored at each time
dir.create("dev/step8", showWarnings = FALSE, recursive = TRUE)
library(ggplot2)
ggsave("dev/step8/synthetic_identifiable_data.png",
       ggplot(d, aes(time, vl, group = id, color = factor(vac))) +
         geom_line(alpha = 0.5) + geom_point(aes(shape = factor(censor))) +
         geom_hline(yintercept = lod, linetype = "dashed") +
         scale_shape_manual(values = c(`0` = 16, `1` = 1), labels = c("Observed", "Censored"), name = NULL) +
         facet_wrap(~ paste0("tlr7=", tlr7, ", vac=", vac)) +
         scale_color_brewer(palette = "Dark2") +
         labs(color = "vac", y = "log10 VL") + theme_report(),
       width = 9, height = 6)

# Check: fit the true model from the true values; are all standard errors estimated? ####
lixoftConnectors::initializeLixoftConnectors(software = "monolix", path = MONOLIX, force = TRUE)
sink("dev/step8/synthetic_identifiable_fim.log", split = TRUE)
create_project(MODEL, file.path(HERE, "synthetic_vl_data.txt"),
               headerTypes = c("id", "time", "observation", "cens", "catcov", "catcov"),
               error_model = list(vl = "constant"), exclude_random_effects = NULL,
               treatment_effects = list(ta = "vac"))
pp <- lixoftConnectors::getPopulationParameterInformation()
truth <- read.csv(file.path(HERE, "true_parameters.csv"))
pp$initialValue[match(truth$parameter, pp$name, nomatch = 0)] <- truth$value[truth$parameter %in% pp$name]
lixoftConnectors::setPopulationParameterInformation(pp)
print(pp)   # initial values = true values (any name not matched keeps Monolix's default)
lixoftConnectors::saveProject("dev/step8/synthetic_identifiable_check.mlxtran")
lixoftConnectors::runScenario()
est <- unlist(lixoftConnectors::getEstimatedPopulationParameters())
se <- lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation
fit <- merge(data.frame(parameter = names(est), estimate = est), truth, by = "parameter", all.x = TRUE)
fit <- merge(fit, se, by = "parameter", all.x = TRUE)
names(fit)[names(fit) == "value"] <- "true"
print(fit, digits = 3)
cat("Condition number:",
    condition_number_of(lixoftConnectors::getCorrelationOfEstimates()$stochasticApproximation), "\n")
# Expect: every rse finite (ideally < 50%) and estimates near the true values
sink()
