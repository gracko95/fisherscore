# Speed comparison: stats::glm() vs fs_glm(engine = "R") vs fs_glm(engine = "cpp")
#
# Run from the package root with:  Rscript inst/benchmarks/benchmark.R
# Writes inst/benchmarks/results.csv and inst/benchmarks/benchmark.png.

library(fisherscore)
stopifnot(requireNamespace("bench", quietly = TRUE))

sim <- function(n, p, family) {
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  beta <- c(-0.5, rep(0.3, p) / sqrt(p))
  eta <- drop(cbind(1, X) %*% beta)
  y <- switch(family,
    binomial = rbinom(n, 1, plogis(eta)),
    poisson  = rpois(n, exp(eta)))
  data.frame(y = y, X)
}

set.seed(651)
grid <- expand.grid(n = c(1e3, 1e4, 1e5, 1e6), p = c(5, 50),
                    family = c("binomial", "poisson"), stringsAsFactors = FALSE)
results <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  d <- sim(g$n, g$p, g$family)
  fam <- get(g$family)()
  iters <- if (g$n >= 1e6) 3 else 10
  b <- bench::mark(
    glm    = coef(glm(y ~ ., fam, d)),
    fs_R   = coef(fs_glm(y ~ ., fam, d, engine = "R")),
    fs_cpp = coef(fs_glm(y ~ ., fam, d, engine = "cpp")),
    iterations = iters, check = function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-8)),
    filter_gc = FALSE
  )
  results[[i]] <- data.frame(g, method = as.character(b$expression),
                             median_sec = as.numeric(b$median),
                             mem_mb = as.numeric(b$mem_alloc) / 2^20)
  message(sprintf("%-8s n=%-7g p=%-3d done", g$family, g$n, g$p))
}
res <- do.call(rbind, results)
out_dir <- file.path("inst", "benchmarks")
write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)
print(res, digits = 3)

png(file.path(out_dir, "benchmark.png"), width = 1000, height = 700, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
cols <- c(glm = "#444444", fs_R = "#1b9e77", fs_cpp = "#d95f02")
for (fam in c("binomial", "poisson")) for (p in c(5, 50)) {
  sub <- res[res$family == fam & res$p == p, ]
  plot(NA, xlim = range(sub$n), ylim = range(sub$median_sec), log = "xy",
       xlab = "n (observations)", ylab = "median time (s)", main = sprintf("%s, p = %d", fam, p))
  for (m in names(cols)) {
    s <- sub[sub$method == m, ]
    lines(s$n, s$median_sec, type = "b", pch = 19, col = cols[m])
  }
  legend("topleft", names(cols), col = cols, lwd = 2, pch = 19, bty = "n")
}
par(op)
dev.off()
