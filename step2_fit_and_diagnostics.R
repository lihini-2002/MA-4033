data_file <- "Sheet1.xlsx"
my_col    <- "220361D"

library(readxl); library(TSA)
dir.create("step2_out", showWarnings = FALSE)
sink("step2_out/step2_output.txt", split = TRUE)

y <- na.omit(ts(as.numeric(read_excel(data_file, sheet = 1)[[my_col]])))
n <- length(y)

summ <- function(fit, name) {
  free <- colnames(fit$var.coef)                          
  k <- length(free) + 1                                   
  if (length(free) == 0) {
    tab <- data.frame(Model = name, Parameter = "(none)", Estimate = NA,
                      SE = NA, t_ratio = NA, Significant = NA)
  } else {
    est <- fit$coef[free]; se <- sqrt(diag(fit$var.coef))[free]
    tab <- data.frame(Model = name, Parameter = free,
                      Estimate = round(est, 4), SE = round(se, 4),
                      t_ratio = round(est / se, 2),
                      Significant = ifelse(abs(est / se) > 2, "yes", "no"))
  }
  crit <- data.frame(Model = name, sigma2 = round(fit$sigma2, 4),
                     logLik = round(fit$loglik, 2), AIC = round(fit$aic, 2),
                     BIC = round(fit$aic - 2 * k + log(n) * k, 2))
  list(tab = tab, crit = crit)
}
root_mod <- function(fit) {          # moduli of AR and MA characteristic roots
  cf <- fit$coef
  ar <- cf[grep("^ar", names(cf))]; ma <- cf[grep("^ma", names(cf))]
  out <- c()
  if (length(ar)) out <- c(out, AR = paste(round(Mod(polyroot(c(1, -ar))), 4), collapse = ", "))
  if (length(ma)) out <- c(out, MA = paste(round(Mod(polyroot(c(1,  ma))), 4), collapse = ", "))
  out
}


fits <- list(
  "ARIMA(0,0,0) with mean" = arima(y, order = c(0,0,0), method = "ML"),
  "ARIMA(0,0,0) zero mean" = arima(y, order = c(0,0,0), method = "ML", include.mean = FALSE),
  "ARIMA(1,0,0)"           = arima(y, order = c(1,0,0), method = "ML"),
  "ARIMA(0,0,1)"           = arima(y, order = c(0,0,1), method = "ML"),
  "ARIMA(1,0,1)"           = arima(y, order = c(1,0,1), method = "ML"),
  "ARIMA(2,0,2)"           = arima(y, order = c(2,0,2), method = "ML"),
  "Subset MA(4): ma4 only" = arima(y, order = c(0,0,4), fixed = c(0,0,0,NA,NA),
                                   transform.pars = FALSE, method = "ML"),
  "ARIMA(0,1,1)"           = arima(y, order = c(0,1,1), method = "ML")
)

all_tab <- NULL; all_crit <- NULL
for (nm in names(fits)) {
  s <- summ(fits[[nm]], nm)
  all_tab <- rbind(all_tab, s$tab); all_crit <- rbind(all_crit, s$crit)
}
cat("==== Table A1: Parameter estimates of candidate models (ML) ====\n")
print(all_tab, row.names = FALSE)
cat("\n==== Table A2: Information criteria of candidate models ====\n")
print(all_crit, row.names = FALSE)
write.csv(all_tab,  "step2_out/tableA1_estimates.csv", row.names = FALSE)
write.csv(all_crit, "step2_out/tableA2_criteria.csv",  row.names = FALSE)

cat("\n==== Table A3: Moduli of characteristic roots (>1 needed) ====\n")
for (nm in c("ARIMA(1,0,0)", "ARIMA(0,0,1)", "ARIMA(1,0,1)", "ARIMA(2,0,2)", "ARIMA(0,1,1)")) {
  cat(nm, ": "); print(root_mod(fits[[nm]]))
}


cat("\n==== Part B: estimates of the mean, CSS vs ML ====\n")
m_ml  <- arima(y, order = c(0,0,0), method = "ML")
m_css <- arima(y, order = c(0,0,0), method = "CSS")
print(m_ml); print(m_css)
cat("Sample mean =", round(mean(y), 4), "  sample variance =", round(var(y), 4), "\n")
cat("\n-- t-test of H0: mu = 0 --\n"); print(t.test(y, mu = 0))


final <- fits[["ARIMA(0,0,0) zero mean"]]
cat("\n==== FINAL MODEL ====\n"); print(final)


res  <- residuals(final)
sres <- res / sqrt(final$sigma2)

png("step2_out/07_std_residuals.png", width = 900, height = 450)
plot(sres, type = "o", pch = 20, cex = 0.6, ylab = "Standardized residuals",
     main = "Standardized residuals - ARIMA(0,0,0)")
abline(h = 0, lty = 1); abline(h = c(-3, 3), lty = 2, col = "red")
dev.off()
cat("\nStd. residuals with |value| > 3:", which(abs(sres) > 3), "\n")

png("step2_out/08_qq_hist.png", width = 900, height = 450)
par(mfrow = c(1, 2))
qqnorm(res, main = "Normal Q-Q plot of residuals"); qqline(res, col = "red")
hist(res, breaks = 20, freq = FALSE, main = "Histogram of residuals", xlab = "Residual")
curve(dnorm(x, 0, sqrt(final$sigma2)), add = TRUE, col = "red", lwd = 2)
dev.off()

cat("\n-- Shapiro-Wilk normality test --\n"); print(shapiro.test(res))
cat("\n-- Runs test (independence) --\n");    print(runs(res))

png("step2_out/09_resid_acf_pacf.png", width = 900, height = 750)
par(mfrow = c(2, 1))
acf(res,  lag.max = 24, main = "Sample ACF of residuals")
pacf(res, lag.max = 24, main = "Sample PACF of residuals")
dev.off()

cat("\n==== Table C1: Ljung-Box tests on residuals ====\n")
npq <- 0                                    
lb <- t(sapply(c(6, 12, 18, 24), function(K) {
  b <- Box.test(res, lag = K, type = "Ljung-Box", fitdf = npq)
  c(K = K, Q_star = round(unname(b$statistic), 3), df = unname(b$parameter),
    p_value = round(b$p.value, 4))
}))
print(lb)
write.csv(lb, "step2_out/tableC1_ljungbox.csv", row.names = FALSE)

png("step2_out/10_tsdiag.png", width = 900, height = 900)
tsdiag(final, gof.lag = 24)
dev.off()


cat("\n==== Part D: overfitting ====\n")
cat("Compare ARIMA(0,0,0) with ARIMA(1,0,0) and ARIMA(0,0,1) in Table A1/A2:\n")
cat("the extra phi / theta should be insignificant (|t| < 2) and AIC should not drop.\n")
print(all_crit[all_crit$Model %in% c("ARIMA(0,0,0) with mean", "ARIMA(1,0,0)", "ARIMA(0,0,1)"), ],
      row.names = FALSE)

cat("\n==== Part E: ARIMA(0,1,1) (what happens if we difference) ====\n")
print(fits[["ARIMA(0,1,1)"]])
cat("Lag-1 sample ACF of diff(y):", round(acf(diff(y), plot = FALSE)$acf[1], 3),
    " (theoretical value -0.5 if y is white noise)\n")

sink()
cat("\nDone. Upload everything in the 'step2_out' folder.\n")
