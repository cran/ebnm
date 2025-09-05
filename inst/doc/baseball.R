## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#", collapse = TRUE, results = "hold",
                      fig.align = "center", dpi = 120)

## ----load-data----------------------------------------------------------------
library("ebnm")
data("wOBA")
nrow(wOBA)
head(wOBA)

## ----woba-histogram, fig.height=2, fig.width=4--------------------------------
library("ggplot2")
ggplot(wOBA, aes(x = x)) +
  geom_histogram(bins = 64, color = "white",fill = "black") +
  theme_classic()

## ----ebnm-normal--------------------------------------------------------------
x <- wOBA$x
s <- wOBA$s
names(x) <- wOBA$Name
names(s) <- wOBA$Name
fit_normal <- ebnm(x, s, prior_family = "normal", mode = "estimate")

## ----ebnm-normal-2------------------------------------------------------------
fit_normal <- ebnm_normal(x, s, mode = "estimate")

## ----summary-ebnm-normal------------------------------------------------------
summary(fit_normal)

## ----plot-ebnm-normal, fig.height=3, fig.width=3.1----------------------------
plot(fit_normal)

## ----plot-ebnm-normal-better, fig.height=3, fig.width=3.7---------------------
plot(fit_normal) +
  geom_point(aes(color = sqrt(wOBA$PA))) +
  labs(x = "wOBA", y = "EB estimate of true wOBA skill", 
       color = expression(sqrt(PA))) +
  scale_color_gradient(low = "blue", high = "red")

## ----fitted-ebnm-normal-------------------------------------------------------
print(head(fitted(fit_normal)), digits = 3)

## ----ebnm-unimodal------------------------------------------------------------
fit_unimodal <- ebnm(x, s, prior_family = "unimodal", mode = "estimate")

## ----ebnm-normal-vs-unimodal, fig.width=4.5, fig.height=2.25------------------
top50 <- order(wOBA$PA, decreasing = TRUE)
top50 <- top50[1:50]
plot(fit_normal, fit_unimodal, subset = top50)

## ----ebnm-normal-vs-unimodal-2------------------------------------------------
dat <- cbind(wOBA[, c("PA","x")],
             fitted(fit_normal),
             fitted(fit_unimodal))
names(dat) <- c("PA", "x", "mean_n", "sd_n", "mean_u", "sd_u")
print(head(dat), digits = 3)

## ----ebnm-normal-vs-unimodal-3, fig.height=2.25, fig.width=6.5, warning=FALSE----
library("cowplot")
p1 <- plot(fit_normal, fit_unimodal, incl_cdf = TRUE, incl_pm = FALSE) +
  xlim(c(.250, .350)) +
  guides(color = "none")
p2 <- plot(fit_normal, fit_unimodal, incl_cdf = TRUE, incl_pm = FALSE) +
  lims(x = c(.350, .450), y = c(0.95, 1))
plot_grid(p1, p2, nrow = 1, ncol = 2, rel_widths = c(1,2))

## ----ebnm-unimodal-confint----------------------------------------------------
fit_unimodal <- ebnm_add_sampler(fit_unimodal)
set.seed(1)
print(head(confint(fit_unimodal, level = 0.8)), digits = 3)

## ----ebnm-npmle---------------------------------------------------------------
fit_npmle <- ebnm(x, s, prior_family = "npmle")

## ----ebnm-npmle-verbose-------------------------------------------------------
fit_npmle <- ebnm(x, s, prior_family = "npmle", 
                  control = list(verbose = TRUE))

## ----plot-npmle, fig.height=2.25, fig.width=4.25------------------------------
plot(fit_normal, fit_unimodal, fit_npmle, incl_cdf = TRUE, subset = top50)

## ----loglik-npmle-------------------------------------------------------------
logLik(fit_unimodal)
logLik(fit_npmle)

## ----ebnm-npmle-finer---------------------------------------------------------
scale_npmle <- ebnm_scale_npmle(x, s, KLdiv_target = 0.001/length(x), 
                                max_K = 1000)
fit_npmle_finer <- ebnm_npmle(x, s, scale = scale_npmle)
logLik(fit_npmle)
logLik(fit_npmle_finer)

## ----ebnm-npmle-sampler-------------------------------------------------------
fit_npmle <- ebnm_add_sampler(fit_npmle)
print(head(quantile(fit_npmle, probs = c(0.1, 0.9))), digits = 3)

## ----ebnm-npmle-confint-------------------------------------------------------
confint(fit_npmle, level = 0.8, parm = "Aaron Judge")

## ----ebnm-deconvolver, fig.height=2.5, fig.width=2.5--------------------------
fit_deconv <- ebnm_deconvolver(x / s, output = ebnm_output_all()) 
plot(fit_deconv, incl_cdf = TRUE, incl_pm = FALSE)

## ----ebnm-deconvolver-head----------------------------------------------------
print(head(confint(fit_deconv, level = 0.8) * s), digits = 3)

## ----session-info-------------------------------------------------------------
sessionInfo()

