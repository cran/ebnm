## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#", collapse = TRUE, results = "hold",
                      fig.align = "center", dpi = 90)

## ----sim-data-----------------------------------------------------------------
set.seed(1)
n <- 400
u <- 2 + (runif(n) < 0.2) * rnorm(n)

## ----sim-data-2---------------------------------------------------------------
s <- rep(1/3, n)
x <- u + s * rnorm(n)

## ----plot-mle, fig.height=3.5, fig.width=3.5----------------------------------
par(mar = c(4, 4, 2, 2))
lims <- c(-0.55, 5.05)
plot(u, x, pch = 4, cex = 0.75, xlim = lims, ylim = lims,
     xlab = "true value", ylab = "estimate", main = "MLE")
abline(a = 0, b = 1, col = "magenta", lty = "dotted")

## ----ebnm-normal--------------------------------------------------------------
library("ebnm")
fit_normal <- ebnm(x, s, prior_family = "normal", mode = "estimate")

## ----plot-ebnm-normal, fig.height=3.5, fig.width=3.5--------------------------
y <- coef(fit_normal)
par(mar = c(4, 4, 2, 2))
plot(u, y, pch = 4, cex = 0.75, xlim = lims, ylim = lims,
     xlab = "true value", ylab = "estimate", main = "normal prior")
abline(a = 0, b = 1, col = "magenta", lty = "dotted")

## ----mse-1--------------------------------------------------------------------
err_mle           <- (x - u)^2
err_shrink_normal <- (y - u)^2
print(round(digits = 4,
            x = c(mle           = sqrt(mean(err_mle)),
                  shrink_normal = sqrt(mean(err_shrink_normal)))))

## ----plot-mse-1, fig.height=3.5, fig.width=3.5--------------------------------
par(mar = c(4, 4, 2, 2))
plot(err_mle, err_shrink_normal, pch = 4, cex = 0.75,
     xlim = c(0, 1.2), ylim = c(0, 1.2))
abline(a = 0, b = 1, col = "magenta", lty = "dotted")

## ----ebnm-pn------------------------------------------------------------------
fit_pn <- ebnm(x, s, prior_family = "point_normal", mode = "estimate")

## ----plot-ebnm-pn, fig.height=3.5, fig.width=3.5------------------------------
par(mar = c(4, 4, 2, 2))
y <- coef(fit_pn)
plot(u, y, pch = 4, cex = 0.75, xlim = lims, ylim = lims,
     xlab = "true value", ylab = "estimate", main = "point-normal prior")
abline(a = 0, b = 1, col = "magenta", lty = "dotted")

## ----mse-2--------------------------------------------------------------------
err_shrink_pn <- (y - u)^2
print(round(digits = 4,
            x = c(mle = sqrt(mean(err_mle)),
                  normal = sqrt(mean(err_shrink_normal)),
                  point_normal = sqrt(mean(err_shrink_pn)))))

## -----------------------------------------------------------------------------
sessionInfo()

