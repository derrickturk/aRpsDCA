\name{bestfit}
\alias{bestfit}
\alias{best.exponential}
\alias{best.hyperbolic}
\alias{best.hyp2exp}
\alias{best.fit}
\alias{best.exponential.from.Np}
\alias{best.hyperbolic.from.Np}
\alias{best.hyp2exp.from.Np}
\alias{best.fit.from.Np}

\title{Best-fitting of Arps decline curves}

\description{
    Perform best-fits of Arps decline curves to rate or cumulative data.
}

\usage{
best.exponential(q, t)
best.hyperbolic(q, t)
best.hyp2exp(q, t)
best.fit(q, t)
best.exponential.from.Np(Np, t)
best.hyperbolic.from.Np(Np, t)
best.hyp2exp.from.Np(Np, t)
best.fit.from.Np(Np, t)
}

\arguments{
  \item{q}{vector of rate data.}
  \item{Np}{vector of cumulative production data.}
  \item{t}{vector of times at which \code{q} or \code{Np} is measured.}
}

\value{
    \code{best.exponential}, \code{best.hyperbolic}, and \code{best.hyp2exp}
    return objects of the appropriate class (as from \code{\link{arps.decline}})
    representing best fits of the appropriate type against \code{q} and
    \code{t}, in the same units as \code{q} and \code{t}.

    \code{best.fit} returns the best overall fit, considering results from each
    function above.

    \code{best.exponential.from.Np}, \code{best.hyperbolic.from.Np}, and
    \code{best.hyp2exp.from.Np} return objects of the appropriate class (as
    from \code{\link{arps.decline}}) representing best fits of the appropriate type
    against \code{Np} and \code{t}, in the same units as \code{Np} and \code{t}.

    \code{best.fit.from.Np} returns the best overall fit, considering results
    from each function above.
}

\details{
    Best-fitting is carried out by minimizing the sum of squared error in the
    rate or cumulative forecast, using \code{\link{nlminb}} as the optimizer.

    Appropriate bounds are applied to decline-curve parameters.
}

\seealso{
  \code{\link{arps}}, \code{\link{nlminb}}
}

\examples{
fitme.hyp2exp.t <- seq(0, 5, 1 / 12) # 5 years
fitme.hyp2exp.q <- hyp2exp.q(
    1000, # Bbl/d
    as.nominal(0.70), # / year
    1.9,
    as.nominal(0.15), # / year
    fitme.hyp2exp.t
) * rnorm(n=length(fitme.hyp2exp.t), mean=1, sd=0.1) # perturb

hyp2exp.fit <- best.hyp2exp(fitme.hyp2exp.q, fitme.hyp2exp.t)
cat(paste("SSE:", hyp2exp.fit$sse))
dev.new()
plot(fitme.hyp2exp.q ~ fitme.hyp2exp.t, main="Hyperbolic-to-Exponential Fit",
     col="blue", log="y", xlab="Time", ylab="Rate")
lines(arps.q(hyp2exp.fit$decline, fitme.hyp2exp.t) ~ fitme.hyp2exp.t,
      col="red")
legend("topright", pch=c(1, NA), lty=c(NA, 1), col=c("blue", "red"), legend=c("Actual", "Fit"))
}
