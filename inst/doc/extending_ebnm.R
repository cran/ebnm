## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#", collapse = TRUE, results = "hold",
                      fig.align = "center", dpi = 90)

## ----tdist--------------------------------------------------------------------
tdist <- function (scale, df) {
  structure(data.frame(scale, df), class = "tdist")
}

## ----opt_t--------------------------------------------------------------------
opt_t <- function (x, s, sigma_init, nu_init) {
  optim(
    par = c(sigma_init, nu_init), 
    fn = function (par) -llik_t(x, s, par[1], par[2]), 
    method = "L-BFGS-B",
    lower = c(min(s)/10, 1),
    upper = c(max(x), 1e3)
  )
}

## ----llik_t-------------------------------------------------------------------
llik_t <- function (x, s, sigma, nu) {
  lik_one_obs <- function (x, s) {
    integrate(lik_times_prior, -Inf, Inf, x = x, s = s,
	          sigma = sigma, nu = nu)$value
  }
  vlik <- Vectorize(lik_one_obs) 
  return(sum(log(vlik(x, s))))
}

## ----lik_times_prior----------------------------------------------------------
lik_times_prior <- function (theta, x, s, sigma, nu) {
   dnorm(x - theta, sd = s) * dt(theta / sigma, df = nu) / sigma
}

## ----opt_grad, eval=FALSE-----------------------------------------------------
# opt_t <- function (x, s, sigma_init, nu_init) {
#   optim(
#     par = c(sigma_init, nu_init),
#     fn = function (par) -llik_t(x, s, par[1], par[2]),
#     gr = function (par) -grad_t(x, s, par[1], par[2]),
#     method = "L-BFGS-B",
#     lower = c(min(s)/10, 1),
#     upper = c(max(x), 1e3)
#   )
# }

## ----grad_t, eval=FALSE-------------------------------------------------------
# library(numDeriv)
# grad_t <- function (x, s, sigma, nu) {
#   grad(function(par) llik_t(x, s, par[1], par[2]), c(sigma, nu))
# }

## ----post_summary-t-----------------------------------------------------------
post_summary_t <- function (x, s, sigma, nu) {
  samp <- post_sampler_t(x, s, sigma, nu, nsamp = 1000)
  return(data.frame(
    mean = colMeans(samp),
    sd = apply(samp, 2, sd),
    second_moment = apply(samp, 2, function (x) mean(x^2))
  ))
}

## ----post_sampler_t-----------------------------------------------------------
# install.packages("mcmc")
library(mcmc)
post_sampler_t <- function (x, s, sigma, nu, nsamp) {
  sample_one_theta <- function (x_i, s_i) {
    lpostdens <- function (theta) {
      dt(theta/sigma, df = nu, log = TRUE) -
	    log(sigma) + 
        dnorm(x_i - theta, sd = s_i, log = TRUE)
    }
    metrop(lpostdens, initial = x_i, nbatch = nsamp)$batch
  }
  vsampler <- Vectorize(sample_one_theta)
  return(vsampler(x, s))
}

## ----ebnm_t-------------------------------------------------------------------
ebnm_t <- function (x, 
                    s = 1, 
					mode = 0, 
					scale = "estimate", 
					g_init = NULL, 
					fix_g = FALSE, 
					output = ebnm_output_default(),
					optmethod = NULL,
					control = NULL) {
				   
  # Some basic argument checks.
  if (mode != 0) {
    stop("The mode of the t-prior must be fixed at zero.")
  }
  if (scale != "estimate") {
    stop("The scale of the t-prior must be estimated rather than fixed ",
	     "at a particular value.")
  }
  
  # If g_init is provided, extract the parameters. Otherwise, provide
  # reasonable initial estimates.
  if (!is.null(g_init)) {
    sigma_init <- g_init$scale
    nu_init    <- g_init$df
  } else {
    sigma_init <- sqrt(mean(x^2))
    nu_init    <- 4
  }
  
  # If g is fixed, use g_init. Otherwise optimize g.
  if (fix_g) {
    sigma <- sigma_init
    nu    <- nu_init
    llik  <- llik_t(x, s, sigma, nu)
  } else {
    opt_res <- opt_t(x, s, sigma_init, nu_init)
    sigma   <- opt_res$par[1]
    nu      <- opt_res$par[2]
    llik    <- -opt_res$value
  }
  
  # Prepare the final output.
  retval <- structure(list(
    data = data.frame(x = x, s = s),
    posterior = post_summary_t(x, s, sigma, nu),
    fitted_g = tdist(scale = sigma, df = nu),
    log_likelihood = llik,
    post_sampler = function (nsamp) post_sampler_t(x, s, sigma, nu, nsamp)
  ), class = c("list", "ebnm"))
  
  return(retval)
}

## ----check--------------------------------------------------------------------
library(ebnm)
set.seed(1)
x <- rnorm(10, sd = 2)
s <- rep(1, 10)
ebnm_check_fn(ebnm_t, x, s)

## ----sim-data-----------------------------------------------------------------
set.seed(1)
theta <- 2 * rt(100, df = 5)
x <- theta + rnorm(100)

## ----t-vs-normal--------------------------------------------------------------
normal_res <- ebnm_normal(x, s = 1)
t_res <- ebnm_t(x, s = 1)

## ----plot-t-vs-normal, fig.width=6, fig.height=4------------------------------
plot(normal_res, t_res)

## ----rmse---------------------------------------------------------------------
rmse_normal <- sqrt(mean((coef(normal_res) - theta)^2))
rmse_t <- sqrt(mean((coef(t_res) - theta)^2))
c(rmse_normal = rmse_normal, rmse_t = rmse_t)

## ----t-fitted-g---------------------------------------------------------------
c(t_res$fitted_g)

## ----session-info-------------------------------------------------------------
sessionInfo()

