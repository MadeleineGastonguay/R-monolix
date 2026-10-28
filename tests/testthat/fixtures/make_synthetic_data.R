# Made-up viral load data with the same columns as the VACTLR7 study data, for tests only.
# Regenerate with: Rscript tests/testthat/fixtures/make_synthetic_data.R
set.seed(42)
ids <- 1:12
times <- c(0, 4, 8, 11, 15, 22, 29, 36)
arms <- data.frame(id = ids,
                   tlr7 = rep(c(0, 1, 0, 1), each = 3),
                   vac = rep(c(0, 0, 1, 1), each = 3),
                   study = rep(c(0, 1), times = 6))
d <- merge(expand.grid(id = ids, time = times), arms)
d <- d[order(d$id, d$time), ]
onset <- 6 + 3 * d$vac + rnorm(nrow(d), 0, 0.5)
raw <- ifelse(d$time < onset, 1, 3 + 2 * (1 - exp(-(d$time - onset) / 5)) - 0.5 * d$vac) + rnorm(nrow(d), 0, 0.3)
lod <- 1.7
d$censor <- as.integer(raw <= lod)
d$vl <- round(ifelse(d$censor == 1, lod, raw), 3)
d$intTLR7VAC <- d$tlr7 * d$vac
d$intTLR7STUD <- d$tlr7 * d$study
d <- d[, c("id", "time", "vl", "censor", "tlr7", "vac", "study", "intTLR7VAC", "intTLR7STUD")]
out <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "synthetic_vl_data.txt")
write.table(d, out, sep = "\t", row.names = FALSE, quote = FALSE)
