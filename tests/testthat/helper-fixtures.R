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
