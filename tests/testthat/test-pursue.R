test_that("abundance arm reproduces limma on complete data without centring or depth adjustment", {
  set.seed(11)
  n <- 30; m <- 40
  meta <- data.frame(group = factor(rep(c("a", "b"), each = n / 2)), age = rnorm(n))
  depth <- rep(1e4, n)
  Y <- matrix(rnorm(n * m, -5, 1), n, m); Y[meta$group == "b", 1:5] <- Y[meta$group == "b", 1:5] + 1
  otu <- round(exp(Y) * depth); otu[otu == 0] <- 1; colnames(otu) <- paste0("f", 1:m)
  design <- PURSUE:::build_design(~ group + age, meta, "group", rowSums(otu))
  arm <- PURSUE:::fit_abundance_arm(otu, rowSums(otu), design, depth_adjust = FALSE, winsorize = FALSE)
  tst <- PURSUE:::abundance_tests(arm, NULL)
  # limma on the same response with the same design
  X <- model.matrix(~ group + age, meta)
  fit <- limma::eBayes(limma::lmFit(t(log(otu / rowSums(otu))), X), robust = TRUE, trend = TRUE)
  tt <- limma::topTable(fit, coef = "groupb", number = Inf, sort.by = "none")
  expect_equal(unname(tst$abund_lfc2 * log(2)), unname(tt$logFC), tolerance = 1e-8)
  expect_equal(unname(tst$abund_p), unname(tt$P.Value), tolerance = 1e-6)
})

test_that("empirical-null centring recovers a known shift", {
  set.seed(2)
  b <- c(rnorm(300, 0.7, 0.1), rnorm(60, 2.5, 0.3)); se <- rep(0.1, 360)
  ce <- PURSUE:::enull_center(b, se)
  expect_lt(abs(ce$delta - 0.7), 0.03)
  expect_gt(ce$pi0, 0.75); expect_lt(ce$pi0, 0.9)
})

test_that("cauchy combination is a valid p-value combiner", {
  P <- cbind(c(0.5, 0.01, NA), c(0.5, 0.5, 0.2))
  cp <- PURSUE:::cauchy_combine(P)
  expect_equal(cp[1], 0.5, tolerance = 1e-8)
  expect_lt(cp[2], 0.05); expect_equal(cp[3], 0.2, tolerance = 1e-8)
})

test_that("pursue runs on sparse data with structural zeros and reports both arms", {
  set.seed(3)
  n <- 40; m <- 80
  x <- rep(0:1, each = n / 2); N <- round(exp(rnorm(n, log(1e4), 0.4)))
  a <- rnorm(m, 0.5, 1); mu <- rnorm(m, log(2e-4), 1.5)
  X <- sapply(1:m, function(j) rbinom(n, 1, plogis(a[j])) * rpois(n, N * exp(mu[j] + rnorm(n, 0, 1))))
  colnames(X) <- paste0("f", 1:m)
  meta <- data.frame(group = factor(x))
  res <- pursue(X, meta, ~ group, "group", min_prevalence = 0, verbose = FALSE)
  r <- res$results
  expect_equal(nrow(r), m)
  expect_true(any(is.finite(r$pres_p)))
  expect_true(any(is.finite(r$abund_p)))
  expect_true(all(r$pres_p[is.finite(r$pres_p)] >= 0 & r$pres_p[is.finite(r$pres_p)] <= 1))
  expect_s3_class(res, "pursue")
})

test_that("input checks fire", {
  otu <- matrix(rpois(40, 5), 8, 5); meta <- data.frame(g = factor(rep(1:2, 4)))
  expect_error(pursue(otu, meta, ~ h, "h", verbose = FALSE), "not in `meta`")
  expect_error(pursue(otu, meta, ~ g, "x", verbose = FALSE), "not a term")
  expect_error(pursue(otu, meta[1:4, , drop = FALSE], ~ g, "g", verbose = FALSE), "one row per sample")
})
