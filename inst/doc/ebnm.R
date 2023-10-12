## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 4, warning = FALSE)

## -----------------------------------------------------------------------------
library(ebnm)
data(wOBA)
nrow(wOBA)
head(wOBA)

## ---- fig.height=3, fig.width=4-----------------------------------------------
library(ggplot2)
ggplot(wOBA, aes(x = x)) +
  geom_histogram(bins = 64, color = "white",fill = "black") +
  theme_classic()

## -----------------------------------------------------------------------------
x <- wOBA$x
s <- wOBA$s
names(x) <- wOBA$Name
names(s) <- wOBA$Name
fit_normal <- ebnm(x, s, prior_family = "normal", mode = "estimate")

## -----------------------------------------------------------------------------
fit_normal <- ebnm_normal(x, s, mode = "estimate")

## -----------------------------------------------------------------------------
summary(fit_normal)

## ---- fig.height=3, fig.width=3.5---------------------------------------------
plot(fit_normal)

## ---- fig.height=3, fig.width=4-----------------------------------------------
plot(fit_normal) +
  geom_point(aes(color = sqrt(wOBA$PA))) +
  labs(x = "wOBA", y = "EB estimate of true wOBA skill", 
       color = expression(sqrt(PA))) +
  scale_color_gradient(low = "blue", high = "red")

## -----------------------------------------------------------------------------
print(head(fitted(fit_normal)), digits = 3)

## -----------------------------------------------------------------------------
fit_unimodal <- ebnm(x, s, prior_family = "unimodal", mode = "estimate")

## ---- fig.width=5.25, fig.height=3--------------------------------------------
top50 <- order(wOBA$PA, decreasing = TRUE)
top50 <- top50[1:50]
plot(fit_normal, fit_unimodal, subset = top50)

## -----------------------------------------------------------------------------
dat <- cbind(wOBA[, c("PA","x")],
             fitted(fit_normal),
             fitted(fit_unimodal))
names(dat) <- c("PA", "x", "mean1", "sd1", "mean2", "sd2")
print(head(dat), digits = 3)

## ---- fig.height=3, fig.width=8-----------------------------------------------
library(cowplot)
p1 <- plot(fit_normal, fit_unimodal, incl_cdf = TRUE, incl_pm = FALSE) +
  xlim(c(.250, .350)) +
  guides(color = "none")
p2 <- plot(fit_normal, fit_unimodal, incl_cdf = TRUE, incl_pm = FALSE) +
  lims(x = c(.350, .450), y = c(0.95, 1))
plot_grid(p1, p2, nrow = 1, ncol = 2, rel_widths = c(3,5))

## -----------------------------------------------------------------------------
fit_unimodal <- ebnm_add_sampler(fit_unimodal)
set.seed(1)
print(head(confint(fit_unimodal, level = 0.8)), digits = 3)

## -----------------------------------------------------------------------------
fit_npmle <- ebnm(x, s, prior_family = "npmle")

## -----------------------------------------------------------------------------
fit_npmle <- ebnm(x, s, prior_family = "npmle", 
                  control = list(verbose = TRUE))

## ---- fig.height=3, fig.width=5-----------------------------------------------
plot(fit_normal, fit_unimodal, fit_npmle, incl_cdf = TRUE, subset = top50)

## -----------------------------------------------------------------------------
logLik(fit_unimodal)
logLik(fit_npmle)

## -----------------------------------------------------------------------------
scale_npmle <- ebnm_scale_npmle(x, s, KLdiv_target = 0.001/length(x), 
                                max_K = 1000)
fit_npmle_finer <- ebnm_npmle(x, s, scale = scale_npmle)
logLik(fit_npmle)
logLik(fit_npmle_finer)

## -----------------------------------------------------------------------------
fit_npmle <- ebnm_add_sampler(fit_npmle)
print(head(quantile(fit_npmle, probs = c(0.1, 0.9))), digits = 3)

## -----------------------------------------------------------------------------
confint(fit_npmle, level = 0.8, parm = "Aaron Judge")

## ---- fig.height=3, fig.width=3-----------------------------------------------
fit_deconv <- ebnm_deconvolver(x / s, output = ebnm_output_all()) 
plot(fit_deconv, incl_cdf = TRUE, incl_pm = FALSE)

## -----------------------------------------------------------------------------
print(head(quantile(fit_deconv, probs = c(0.1, 0.9)) * s), digits = 3)

## -----------------------------------------------------------------------------
sessionInfo()

