# Small fake IGS output: 4 runs x 2 parameters, same shape as initial_guess_sampling()$tabestimates
fake_tabestimates <- function() {
  data.frame(
    parameter = rep(c("k_pop", "omega_k"), times = 4),
    estimate  = c(1.0, 0.5,  1.2, 0.6,  3.0, 0.9,  10, 2),
    rses      = c(10, 20, 12, 22, 40, 50, 200, 300),
    run       = rep(1:4, each = 2),
    bicc      = rep(c(100, 100.5, 104, 130), each = 2)
  )
}

# A profile_likelihood() table: two estimated parameters per grid value
fake_profile <- function(parameter, values, ofv) {
  data.frame(parameter = rep(c("x_pop", "omega_x"), length(values)), estimate = 1, profiled_parameter = parameter,
             profiled_value = rep(values, each = 2), ofv = rep(ofv, each = 2))
}

# A made-up SCM: `models` is a list of covariate models (parameter = covariates), one per entry
fake_mb <- function(models, tested, bicc, best) {
  prm <- c("ta", "NPs", "m"); covs <- c("tlr7", "vac")
  Map(function(m, t, b, i) {
    im <- matrix(FALSE, length(prm), length(covs), dimnames = list(prm, covs))
    for (p in names(m)) im[p, m[[p]]] <- TRUE
    list(tested = t, BICc = b, LL = -b, individualModels = as.data.frame(im), bestModel = i == best)
  }, models, tested, bicc, seq_along(models))
}

# forward: ta-vac (-10), m-tlr7 (-5), NPs-tlr7 (-1), then NPs-vac doesn't help; backward from
# the 3-covariate model removes ta-vac. scm_delta2 stops after m-tlr7, so it keeps ta-vac.
fake_mb_backward <- function() {
  fake_mb(
    models = list(list(), list(ta = "vac"), list(NPs = "tlr7"),
                  list(ta = "vac", m = "tlr7"), list(ta = "vac", NPs = "vac"),
                  list(ta = "vac", m = "tlr7", NPs = "tlr7"),
                  list(ta = "vac", m = "tlr7", NPs = c("tlr7", "vac")),
                  list(m = "tlr7", NPs = "tlr7"), list(ta = "vac", NPs = "tlr7")),
    tested = list(c("", ""), c("ta", "vac"), c("NPs", "tlr7"), c("m", "tlr7"), c("NPs", "vac"),
                  c("NPs", "tlr7"), c("NPs", "vac"), c("ta", "vac"), c("m", "tlr7")),
    bicc = c(100, 90, 95, 85, 89, 84, 86, 83, 88), best = 8)
}

# A candidate's IGS results for compare_models(): BICc per run, optional treatment effects,
# Wald p-values and condition numbers
fake_candidate <- function(bicc, betas = character(), p = NULL, cn = NULL, runs = 1:4) {
  te <- rbind(
    data.frame(parameter = "k_pop", estimate = 1, run = runs),
    data.frame(parameter = "omega_k", estimate = 0.3 + bicc / 1000, run = runs),
    if (length(betas)) data.frame(parameter = rep(betas, each = length(runs)), estimate = 0.1, run = runs)
  )
  te$rse <- 10
  te$bicc <- bicc[match(te$run, runs)]
  out <- list(tabestimates = te)
  if (!is.null(p)) out$pvals <- data.frame(parameter = betas[1], p.value = p, run = runs)
  if (!is.null(cn)) out$condition_number <- data.frame(run = runs, condition_number = cn)
  out
}
