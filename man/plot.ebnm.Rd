% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebnm_methods.R
\name{plot.ebnm}
\alias{plot.ebnm}
\title{Plot an ebnm object}
\usage{
\method{plot}{ebnm}(x, remove_abline = FALSE, ...)
}
\arguments{
\item{x}{The fitted \code{ebnm} object.}

\item{remove_abline}{To better illustrate shrinkage effects, the plot
will include the line \eqn{y = x} by default. If
\code{remove_abline = TRUE}, then this line will not be drawn.}

\item{...}{Additional parameters to be passed to \code{ggplot2} function
\code{\link[ggplot2]{geom_point}}.}
}
\value{
A \code{ggplot} object.
}
\description{
Given a fitted \code{ebnm} object, produces a plot of posterior means vs.
  observations.
}
\details{
An object of class \code{ggplot} is returned, so that the plot can be
  customized in the usual \code{\link[ggplot2]{ggplot2}} fashion.
}
\examples{
theta <- c(rep(0, 100), rexp(100))
s <- 1
x <- theta + rnorm(200, 0, s)
ebnm.res <- ebnm(x, s)
plot(ebnm.res)

# Customize plot:
library(ggplot2)
plot(ebnm.res, color = "blue", remove_abline = TRUE) +
  theme_bw() +
  labs(x = "Simulated data")

}
