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
\alias{best.exponential.from.interval}
\alias{best.hyperbolic.from.interval}
\alias{best.hyp2exp.from.interval}
\alias{best.fit.from.interval}
\alias{best.exponential.curtailed}
\alias{best.hyperbolic.curtailed}
\alias{best.hyp2exp.curtailed}
\alias{best.curtailed.fit}
\alias{best.exponential.curtailed.from.Np}
\alias{best.hyperbolic.curtailed.from.Np}
\alias{best.hyp2exp.curtailed.from.Np}
\alias{best.curtailed.fit.from.Np}
\alias{best.exponential.curtailed.from.interval}
\alias{best.hyperbolic.curtailed.from.interval}
\alias{best.hyp2exp.curtailed.from.interval}
\alias{best.curtailed.fit.from.interval}
\alias{best.exponential.with.buildup}
\alias{best.hyperbolic.with.buildup}
\alias{best.hyp2exp.with.buildup}
\alias{best.fit.with.buildup}
\alias{best.exponential.from.Np.with.buildup}
\alias{best.hyperbolic.from.Np.with.buildup}
\alias{best.hyp2exp.from.Np.with.buildup}
\alias{best.fit.from.Np.with.buildup}
\alias{best.exponential.from.interval.with.buildup}
\alias{best.hyperbolic.from.interval.with.buildup}
\alias{best.hyp2exp.from.interval.with.buildup}
\alias{best.fit.from.interval.with.buildup}

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
    max(q) * 5, # qi < qmax * 5
    10) # = 0.99995 / [time] effective
  )

best.hyperbolic(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
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
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35) # Df <= 0.35
  )

best.exponential.curtailed(q, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
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
    max(q) * 5, # qi < qmax * 5
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
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )

best.fit(q, t)

best.curtailed.fit(q, t)

best.exponential.from.Np(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10) # = 0.99995 / [time] effective)
  )

best.exponential.from.interval(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10) # = 0.99995 / [time] effective)
  )

best.hyperbolic.from.Np(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2)  # b <= 2.0
  )

best.hyperbolic.from.interval(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
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
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35) # Df <= 0.35
  )

best.hyp2exp.from.interval(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35) # Df <= 0.35
  )

best.exponential.curtailed.from.Np(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    t[length(t)])
  )

best.exponential.curtailed.from.interval(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
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
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    t[length(t)])
  )

best.hyperbolic.curtailed.from.interval(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
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
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )

best.hyp2exp.curtailed.from.interval(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0
  ),
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )

best.fit.from.Np(Np, t)

best.fit.from.interval(volume, t, t.begin=0.0)

best.curtailed.fit.from.Np(Np, t)

best.curtailed.fit.from.interval(volume, t, t.begin=0.0)

best.exponential.with.buildup(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10),        # = 0.99995 / [time] effective
  initial.rate=q[1], time.to.peak=t[which.max(q)])

best.hyperbolic.with.buildup(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2), # b <= 2.0
  initial.rate=q[1], time.to.peak=t[which.max(q)])

best.hyp2exp.with.buildup(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=q[1], time.to.peak=t[which.max(q)])

best.fit.with.buildup(q, t)

best.exponential.from.Np.with.buildup(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10), # = 0.99995 / [time] effective
  initial.rate=Np[1] / t[1],
  time.to.peak=(t[which.max(diff(Np))] + t[which.max(diff(Np)) + 1]) / 2.0)

best.exponential.from.interval.with.buildup(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10), # = 0.99995 / [time] effective
  initial.rate=volume[1] / (t[1] - t.begin),
  time.to.peak=(t - diff(c(t.begin, t)) / 2)[which.max(volume)])

best.hyperbolic.from.Np.with.buildup(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2), # b <= 2.0
  initial.rate=Np[1] / t[1],
  time.to.peak=(t[which.max(diff(Np))] + t[which.max(diff(Np)) + 1]) / 2.0)

best.hyperbolic.from.interval.with.buildup(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2), # b <= 2.0
  initial.rate=volume[1] / (t[1] - t.begin),
  time.to.peak=(t - diff(c(t.begin, t)) / 2)[which.max(volume)])

best.hyp2exp.from.Np.with.buildup(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=Np[1] / t[1],
  time.to.peak=(t[which.max(diff(Np))] + t[which.max(diff(Np)) + 1]) / 2.0)

best.hyp2exp.from.interval.with.buildup(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=volume[1] / (t[1] - t.begin),
  time.to.peak=(t - diff(c(t.begin, t)) / 2)[which.max(volume)])

best.fit.from.Np.with.buildup(Np, t)

best.fit.from.interval.with.buildup(volume, t, t.begin=0.0)
}

\arguments{
  \item{q}{vector of rate data.}
  \item{Np}{vector of cumulative production data.}
  \item{volume}{vector of interval volume data.}
  \item{t}{vector of times at which \code{q}, \code{Np}, or \code{volume} is measured.}
  \item{t.begin}{initial time for interval volume data, if non-zero.}
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

    \code{best.exponential.from.interval}, \code{best.hyperbolic.from.interval},
    and \code{best.hyp2exp.from.interval} return objects of the appropriate
    class (as from \code{\link{arps.decline}}) representing best fits of the
    appropriate type against \code{volume} and \code{t}, in the same units as
    \code{volume} and \code{t}.

    For these functions, \code{t} is taken to represent the time at the end of
    each producing interval; the beginning time for the first interval may be
    specified as \code{t.begin} if it is non-zero.

    \code{best.fit.from.interval} returns the best overall fit, considering
    results from each function above.

    \code{best.exponential.curtailed}, \code{best.hyperbolic.curtailed},
    \code{best.hyp2exp.curtailed}, \code{best.curtailed.fit},
    \code{best.exponential.curtailed.from.Np},
    \code{best.hyperbolic.curtailed.from.Np},
    \code{best.hyp2exp.curtailed.from.Np}, \code{best.curtailed.fit.from.Np},
    \code{best.exponential.curtailed.from.interval},
    \code{best.hyperbolic.curtailed.from.interval},
    \code{best.hyp2exp.curtailed.from.interval}, and
    \code{best.curtailed.fit.from.interval} work as the corresponding functions
    above, but may return curtailed declines (as from \code{\link{curtail}}).

    \code{best.exponential.with.buildup}, \code{best.hyperbolic.with.buildup},
    \code{best.hyp2exp.with.buildup}, \code{best.fit.with.buildup},
    \code{best.exponential.from.Np.with.buildup},
    \code{best.hyperbolic.from.Np.with.buildup},
    \code{best.hyp2exp.from.Np.with.buildup},
    \code{best.fit.from.Np.with.buildup},
    \code{best.exponential.from.interval.with.buildup},
    \code{best.hyperbolic.from.interval.with.buildup},
    \code{best.hyp2exp.from.interval.with.buildup}, and
    \code{best.fit.from.interval.with.buildup} work as the corresponding
    functions above, but will return a fit including a linear buildup
    portion (as from \code{\link{arps.with.buildup}}).
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
