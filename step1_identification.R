data_file <- "Sheet1.xlsx"     
my_col    <- "220361D"         


pkgs <- c("readxl", "TSA", "tseries", "fUnitRoots")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(readxl); library(TSA); library(tseries); library(fUnitRoots)

dir.create("step1_out", showWarnings = FALSE)
sink("step1_out/step1_output.txt", split = TRUE)   


dat <- read_excel(data_file, sheet = 1)
y   <- ts(as.numeric(dat[[my_col]]))
y   <- na.omit(y)
cat("Column used:", my_col, "\n")
cat("n =", length(y), "\n\n")
print(summary(as.numeric(y)))
cat("SD =", sd(y), "\n\n")

png("step1_out/01_timeplot.png", width = 900, height = 450)
plot(y, type = "o", pch = 20, cex = 0.6, ylab = "Y_t", xlab = "Time",
     main = "Time series plot of Dataset 1")
abline(h = mean(y), lty = 2, col = "red")
dev.off()

png("step1_out/02_acf_pacf_original.png", width = 900, height = 750)
par(mfrow = c(2, 1))
acf(y,  lag.max = 24, main = "Sample ACF of Dataset 1")    
pacf(y, lag.max = 24, main = "Sample PACF of Dataset 1")
dev.off()

cat("==== Sample ACF (original) ====\n")
print(round(acf(y,  lag.max = 24, plot = FALSE)$acf[, , 1], 3))
cat("\n==== Sample PACF (original) ====\n")
print(round(pacf(y, lag.max = 24, plot = FALSE)$acf[, , 1], 3))
cat("Approx. 2 s.e. bound = +/-", round(2 / sqrt(length(y)), 3), "\n\n")

cat("==== Sample EACF (original) ====\n")
eacf(y, ar.max = 7, ma.max = 13)


k <- ar(diff(y), order.max = 10)$order
cat("\n==== ADF test ====\nk selected by AIC on diff(y):", k, "\n")
cat("\n-- H1: stationary with unknown mean (type = 'c') --\n")
print(adfTest(y, lags = k, type = "c"))
cat("\n-- H1: stationary around a linear trend (type = 'ct') --\n")
print(adfTest(y, lags = k, type = "ct"))
cat("\n-- tseries::adf.test (trend alternative, default lag) --\n")
print(adf.test(y))


dy <- diff(y)
png("step1_out/03_timeplot_diff.png", width = 900, height = 450)
plot(dy, type = "o", pch = 20, cex = 0.6, ylab = expression(nabla*Y[t]),
     main = "First difference of Dataset 1")
abline(h = 0, lty = 2, col = "red")
dev.off()

png("step1_out/04_acf_pacf_diff.png", width = 900, height = 750)
par(mfrow = c(2, 1))
acf(dy,  lag.max = 24, main = "Sample ACF of first difference")
pacf(dy, lag.max = 24, main = "Sample PACF of first difference")
dev.off()

cat("\n==== Sample EACF (first difference) ====\n")
eacf(dy, ar.max = 7, ma.max = 13)

cat("\n==== ADF test on first difference ====\n")
k2 <- ar(diff(dy), order.max = 10)$order
print(adfTest(dy, lags = k2, type = "c"))


png("step1_out/05_armasubsets_original.png", width = 900, height = 600)
plot(armasubsets(y = y, nar = 7, nma = 7, y.name = "Y", ar.method = "ols"))
title("Best subset ARMA (BIC) - original series", line = 3)
dev.off()

png("step1_out/06_armasubsets_diff.png", width = 900, height = 600)
plot(armasubsets(y = dy, nar = 7, nma = 7, y.name = "dY", ar.method = "ols"))
title("Best subset ARMA (BIC) - first difference", line = 3)
dev.off()


grid <- expand.grid(p = 0:3, d = 0:1, q = 0:3)
grid$AIC <- grid$BIC <- grid$AICc <- NA
for (i in seq_len(nrow(grid))) {
  fit <- tryCatch(arima(y, order = c(grid$p[i], grid$d[i], grid$q[i]),
                        method = "ML", include.mean = (grid$d[i] == 0)),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    npar <- length(fit$coef) + 1
    n    <- length(y) - grid$d[i]
    grid$AIC[i]  <- fit$aic
    grid$BIC[i]  <- fit$aic - 2 * npar + log(n) * npar
    grid$AICc[i] <- fit$aic + 2 * npar * (npar + 1) / (n - npar - 1)
  }
}
cat("\n==== Candidate ARIMA models sorted by AIC ====\n")
print(head(grid[order(grid$AIC), ], 12), row.names = FALSE, digits = 6)
cat("\n==== Candidate ARIMA models sorted by BIC ====\n")
print(head(grid[order(grid$BIC), ], 12), row.names = FALSE, digits = 6)

sink()
cat("\nDone. Upload everything in the 'step1_out' folder.\n")
