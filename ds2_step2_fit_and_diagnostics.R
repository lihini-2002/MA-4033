data_file <- "sheet2.xlsx"
my_col    <- "220361D"

library(readxl); library(TSA)
dir.create("ds2_step2_out", showWarnings = FALSE)
sink("ds2_step2_out/ds2_step2_output.txt", split = TRUE)

y  <- na.omit(ts(as.numeric(read_excel(data_file, sheet = 1)[[my_col]])))
dy <- diff(y)
n  <- length(dy)                       

summ <- function(fit, name, pq) {
  free <- colnames(fit$var.coef)
  k <- length(free) + 1
  est <- fit$coef[free]; se <- sqrt(diag(fit$var.coef))[free]
  tab <- data.frame(Model = name, Parameter = free,
                    Estimate = round(est, 4), SE = round(se, 4),
                    t_ratio = round(est / se, 2),
                    Significant = ifelse(abs(est / se) > 2, "yes", "no"))
  r <- residuals(fit)
  lbp <- sapply(c(6, 12, 18, 24), function(K)
    round(Box.test(r, lag = K, type = "Ljung-Box", fitdf = pq)$p.value, 3))
  crit <- data.frame(Model = name, sigma2 = round(fit$sigma2, 4),
                     logLik = round(fit$loglik, 2), AIC = round(fit$aic, 2),
                     AICc = round(fit$aic + 2 * k * (k + 1) / (n - k - 1), 2),
                     BIC = round(fit$aic - 2 * k + log(n) * k, 2),
                     LB6 = lbp[1], LB12 = lbp[2], LB18 = lbp[3], LB24 = lbp[4])
  list(tab = tab, crit = crit)
}
root_mod <- function(fit) {
  cf <- fit$coef
  ar <- cf[grep("^ar", names(cf))]; ma <- cf[grep("^ma", names(cf))]
  out <- c()
  if (length(ar)) out <- c(out, AR = paste(round(Mod(polyroot(c(1, -ar))), 4), collapse = ", "))
  if (length(ma)) out <- c(out, MA = paste(round(Mod(polyroot(c(1,  ma))), 4), collapse = ", "))
  out
}

spec <- list(
  "ARIMA(0,1,0)" = c(0,1,0), "ARIMA(0,1,1)" = c(0,1,1), "ARIMA(1,1,0)" = c(1,1,0),
  "ARIMA(2,1,0)" = c(2,1,0), "ARIMA(0,1,2)" = c(0,1,2), "ARIMA(1,1,1)" = c(1,1,1),
  "ARIMA(1,1,2)" = c(1,1,2), "ARIMA(2,1,1)" = c(2,1,1), "ARIMA(3,1,0)" = c(3,1,0))
fits <- lapply(spec, function(o) arima(y, order = o, method = "ML"))

fits[["Subset ARIMA(2,1,0), phi1=0"]] <- arima(y, order = c(2,1,0), fixed = c(0, NA),
                                               transform.pars = FALSE, method = "ML")
pq <- c(sapply(spec, function(o) o[1] + o[3]), 1)

all_tab <- NULL; all_crit <- NULL
for (i in seq_along(fits)) {
  nm <- names(fits)[i]
  if (length(fits[[nm]]$coef) == 0) {   # random walk
    r <- residuals(fits[[nm]]); k <- 1
    lbp <- sapply(c(6,12,18,24), function(K) round(Box.test(r, K, "Ljung-Box")$p.value, 3))
    all_tab  <- rbind(all_tab, data.frame(Model = nm, Parameter = "(none)", Estimate = NA,
                                          SE = NA, t_ratio = NA, Significant = NA))
    all_crit <- rbind(all_crit, data.frame(Model = nm, sigma2 = round(fits[[nm]]$sigma2, 4),
                      logLik = round(fits[[nm]]$loglik, 2), AIC = round(fits[[nm]]$aic, 2),
                      AICc = round(fits[[nm]]$aic + 2*k*(k+1)/(n-k-1), 2),
                      BIC = round(fits[[nm]]$aic - 2*k + log(n)*k, 2),
                      LB6 = lbp[1], LB12 = lbp[2], LB18 = lbp[3], LB24 = lbp[4]))
  } else {
    s <- summ(fits[[nm]], nm, pq[i])
    all_tab <- rbind(all_tab, s$tab); all_crit <- rbind(all_crit, s$crit)
  }
}
cat("==== Table A1: Parameter estimates (ML) ====\n");   print(all_tab,  row.names = FALSE)
cat("\n==== Table A2: Criteria + Ljung-Box p-values ====\n"); print(all_crit[order(all_crit$AIC), ], row.names = FALSE)
write.csv(all_tab,  "ds2_step2_out/tableA1_estimates.csv", row.names = FALSE)
write.csv(all_crit, "ds2_step2_out/tableA2_criteria.csv",  row.names = FALSE)

cat("\n==== Table A3: Moduli of characteristic roots (>1 needed) ====\n")
for (nm in names(fits)[-1]) { cat(nm, ": "); print(root_mod(fits[[nm]])) }


cat("\n==== Part B1: is a drift (constant in the differenced model) needed? ====\n")
print(t.test(dy, mu = 0))
drift <- arima(dy, order = c(2,0,0), method = "ML")        # ARIMA(2,1,0) + drift
cat("\nARIMA(2,1,0) with drift (fitted to the differences):\n"); print(drift)
cat("t-ratio of drift:", round(drift$coef["intercept"] / sqrt(drift$var.coef["intercept","intercept"]), 2), "\n")

final_order <- c(2,1,0)
final <- fits[["ARIMA(2,1,0)"]]
cat("\n==== Part B2: ARIMA(2,1,0) - ML vs CSS ====\n")
print(final)
print(arima(y, order = final_order, method = "CSS"))

cat("\n==== FINAL MODEL ====\n"); print(final)

res  <- residuals(final)[-1]                 
sres <- res / sqrt(final$sigma2)

png("ds2_step2_out/07_std_residuals.png", width = 900, height = 450)
plot(ts(sres, start = 2), type = "o", pch = 20, cex = 0.6, ylab = "Standardized residuals",
     xlab = "Time", main = "Standardized residuals - ARIMA(2,1,0)")
abline(h = 0); abline(h = c(-3, 3), lty = 2, col = "red")
dev.off()
cat("\nStd. residuals with |value| > 3 at times:", which(abs(sres) > 3) + 1, "\n")
cat("Range of std. residuals:", round(range(sres), 2), "\n")

png("ds2_step2_out/08_qq_hist.png", width = 900, height = 450)
par(mfrow = c(1, 2))
qqnorm(res, main = "Normal Q-Q plot of residuals"); qqline(res, col = "red")
hist(res, breaks = 20, freq = FALSE, main = "Histogram of residuals", xlab = "Residual")
curve(dnorm(x, 0, sqrt(final$sigma2)), add = TRUE, col = "red", lwd = 2)
dev.off()

cat("\n-- Shapiro-Wilk normality test --\n"); print(shapiro.test(res))
cat("\n-- Runs test (independence) --\n");    print(runs(res))

png("ds2_step2_out/09_resid_acf_pacf.png", width = 900, height = 750)
par(mfrow = c(2, 1))
acf(res,  lag.max = 24, main = "Sample ACF of residuals - ARIMA(2,1,0)")
pacf(res, lag.max = 24, main = "Sample PACF of residuals - ARIMA(2,1,0)")
dev.off()

cat("\n==== Table C1: Ljung-Box tests on residuals (df = K - p - q) ====\n")
lb <- t(sapply(c(6, 12, 18, 24), function(K) {
  b <- Box.test(res, lag = K, type = "Ljung-Box", fitdf = 2)
  c(K = K, Q_star = round(unname(b$statistic), 3), df = unname(b$parameter),
    p_value = round(b$p.value, 4))
}))
print(lb)
write.csv(lb, "ds2_step2_out/tableC1_ljungbox.csv", row.names = FALSE)

png("ds2_step2_out/10_tsdiag.png", width = 900, height = 900)
tsdiag(final, gof.lag = 24)
dev.off()

cat("\n==== Part D: overfitting ARIMA(2,1,0) ====\n")
for (nm in c("ARIMA(2,1,0)", "ARIMA(3,1,0)", "ARIMA(2,1,1)")) {
  cat("\n--", nm, "--\n"); print(fits[[nm]])
}
print(all_crit[all_crit$Model %in% c("ARIMA(2,1,0)", "ARIMA(3,1,0)", "ARIMA(2,1,1)"), ],
      row.names = FALSE)

sink()
cat("\nDone. Upload everything in the 'ds2_step2_out' folder.\n")
