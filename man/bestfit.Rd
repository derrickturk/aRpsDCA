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
\alias{best.exponential.curtailed}
\alias{best.hyperbolic.curtailed}
\alias{best.hyp2exp.curtailed}
\alias{best.curtailed.fit}
\alias{best.exponential.curtailed.from.Np}
\alias{best.hyperbolic.curtailed.from.Np}
\alias{best.hyp2exp.curtailed.from.Np}
\alias{best.curtailed.fit.from.Np}
\alias{best.exponential.with.buildup}
\alias{best.hyperbolic.with.buildup}
\alias{best.hyp2exp.with.buildup}
\alias{best.fit.with.buildup}

\title{Best-fitting of Arps decline curves}

\description{
    Perform best-fits of Arps decline curves to rate or cumulative data.
}

\usage{
best.exponential(q, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10) # = 0.99995 / [time] effective
  )

best.hyperbolic(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    2)  # b <= 2.0
  )

best.hyp2exp(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35) # Df <= 0.35
  )

best.fit(q, t)

best.exponential.from.Np(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(Np) * 3, # qi < max(Np) * 3
    10) # = 0.99995 / [time] effective)
  )

best.hyperbolic.from.Np(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(Np) * 3, # qi < max(Np) * 3
    10, # = 0.99995 / [time] effective
    2)  # b <= 2.0
  )

best.hyp2exp.from.Np(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(Np) * 3, # qi < max(Np) * 3
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35) # Df <= 0.35
  )

best.fit.from.Np(Np, t)

best.exponential.curtailed(q, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    t[length(t)])
  )

best.hyperbolic.curtailed(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    t[length(t)])
  )

best.hyp2exp.curtailed(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )

best.curtailed.fit(q, t)

best.exponential.curtailed.from.Np(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(Np) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    t[length(t)])
  )

best.hyperbolic.curtailed.from.Np(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(Np) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    t[length(t)])
  )

best.hyp2exp.curtailed.from.Np(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0
  ),
  upper=c( # upper bounds
    max(Np) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )

best.curtailed.fit.from.Np(Np, t)

best.exponential.with.buildup(q, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0, # initial rate
    0), # time to peak > 0
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10), # = 0.99995 / [time] effective
  initial.rate=q[1], time.to.peak=t[which.max(q)])

best.hyperbolic.with.buildup(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0,  # initial rate
    0), # time to peak > 0
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    2),  # b <= 2.0
  initial.rate=q[1], time.to.peak=t[which.max(q)])

best.hyp2exp.with.buildup(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0,  # initial rate
    0), # time to peak > 0
  upper=c( # upper bounds
    max(q) * 3, # qi < qmax * 3
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=q[1], time.to.peak=t[which.max(q)])

best.fit.with.buildup(q, t)
}

\arguments{
  \item{q}{vector of rate data.}
  \item{Np}{vector of cumulative production data.}
  \item{t}{vector of times at which \code{q} or \code{Np} is measured.}
  \item{lower}{lower bounds for decline parameters (sane defaults are provided).}
  \item{upper}{upper bounds for decline parameters (sane defaults are provided).}
  \item{initial.rate}{initial rate, for declines with buildup.}
  \item{time.to.peak}{time to peak rate, for declines with buildup.}
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

    \code{best.exponential.curtailed}, \code{best.hyperbolic.curtailed},
    \code{best.hyp2exp.curtailed}, \code{best.curtailed.fit},
    \code{best.exponential.curtailed.from.Np},
    \code{best.hyperbolic.curtailed.from.Np},
    \code{best.hyp2exp.curtailed.from.Np}, and \code{best.curtailed.fit.from.Np}
    work as the corresponding functions above, but may return curtailed
    declines (as from \code{\link{curtail}}).

    \code{best.exponential.with.buildup}, \code{best.hyperbolic.with.buildup},
    \code{best.hyp2exp.with.buildup}, and \code{best.fit.with.buildup} work as
    the corresponding functions above, but will return a fit including a linear
    buildup portion (as from \code{\link{arps.with.buildup}}).
}

\details{
    Best-fitting is carried out by minimizing the sum of squared error in the
    rate or cumulative forecast, using \code{\link{nlminb}} as the optimizer.

    Appropriate bounds are applied to decline-curve parameters by default, but
    may be altered using the \code{lower} and \code{upper} arguments to each
    specific function.
}

\seealso{
  \code{\link{arps}}, \code{\link{curtailed}}, \code{\link{arps.with.buildup}},
  \code{\link{nlminb}}
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
