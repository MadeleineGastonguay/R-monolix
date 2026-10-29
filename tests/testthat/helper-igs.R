# Fake initial_guess_sampling() output with every record: 8 runs, a one-compartment-like model
# with k_pop, v_pop, omega_k, omega_v and error parameter a; 3 individuals.
fake_igs_full <- function(n = 8) {
  set.seed(42)
  params <- c("k_pop", "v_pop", "omega_k", "omega_v", "a")
  iv <- do.call(rbind, lapply(seq_len(n), function(r) data.frame(
    parameter = params, initialValue_dummy = NA, run = r,
    initial_value = c(10^runif(1, -1, 1), 10^runif(1, 0, 2), runif(2, 1, 5), 1), method = "MLE")))
  iv$initialValue_dummy <- NULL
  te <- do.call(rbind, lapply(seq_len(n), function(r) data.frame(
    parameter = params, estimate = c(0.5, 20, 0.3, 0.2, 0.4) * exp(rnorm(5, 0, 0.05)),
    rse = c(10, 15, 30, 40, 8), run = r)))
  bic <- data.frame(run = seq_len(n), bicc = 100 + sort(rexp(n, 0.5)))
  te <- merge(te, bic)
  te$ofv <- te$bicc - 10; te$aic <- te$bicc - 5; te$bic <- te$bicc - 1
  tabiters <- do.call(rbind, lapply(seq_len(n), function(r) {
    it <- data.frame(k_pop = seq(1, 0.5, length.out = 20), v_pop = seq(10, 20, length.out = 20),
                     omega_k = 1, omega_v = 1, a = 1, convergenceIndicator = seq(-50, -40, length.out = 20))
    it$iteration <- 1:20; it$run <- r; it
  }))
  phase <- do.call(rbind, lapply(seq_len(n), function(r) c(10, 20)))
  fit <- expand.grid(id = 1:3, time = c(0, 1, 2, 4, 8), run = seq_len(n))
  fit$y <- 20 * exp(-0.5 * fit$time) * (1 + fit$id / 10)
  fit$popPred <- 20 * exp(-0.5 * fit$time)
  fit$indivPred_mode <- fit$y * 1.02
  list(tabestimates = te, initialvalues = iv, tabiters = tabiters, phase = phase, indivFit = fit,
       condition_number = data.frame(run = seq_len(n), condition_number = 50),
       shrinkage = data.frame(run = rep(seq_len(n), each = 2), parameter = c("k", "v"),
                              conditionalMean = 10, conditionalMode = 12, conditionalDist = 10))
}
