\name{hyperbolic}
\alias{hyperbolic}
\alias{hyperbolic.q}
\alias{hyperbolic.Np}
\alias{hyperbolic.D}

\title{Arps hyperbolic declines}

\description{
    Compute rates, cumulative production values, and instantaneous nominal declines for Arps hyperbolic decline curves.
}

\usage{
hyperbolic.q(qi, Di, b, t)
hyperbolic.Np(qi, Di, b, t)
hyperbolic.D(Di, b, t)
}

\arguments{
  \item{qi}{initial rate [volume / time], i.e. q(t = 0).}
  \item{Di}{initial nominal Arps decline exponent [1 / time].}
  \item{b}{Arps hyperbolic exponent.}
  \item{t}{time at which to evaluate rate or cumulative [time].}
}

\value{
    \code{hyperbolic.q} returns the rate for each element of \code{t},
    in the same units as \code{qi}.

    \code{hyperbolic.Np} returns the cumulative production for each element of
    \code{t}, in the same units as \code{qi * t}.

    \code{hyperbolic.D} returns the nominal instantaneous decline for each
    element of \code{t}. This can be converted to effective decline and
    rescaled in time by use of \code{\link{as.effective}} and
    \code{\link{rescale.by.time}}.
}

\details{
    Assumes consistent units of time between \code{qi}, \code{D}, and \code{t}. To convert, see the decline-rate conversion functions referenced below.

    When appropriate, internally uses \code{\link{harmonic.q}} and \code{\link{harmonic.Np}} to avoid singularities in calculations for \code{b} near 1.
}

\seealso{
  \code{\link{as.effective}}, \code{\link{as.nominal}}, \code{\link{rescale.by.time}}.
}

\examples{
## qi = 1000 Mscf/d, Di = 95% effective / year, b = 1.2, t from 0 to 25 days
hyperbolic.q(1000, as.nominal(0.95, from.period="year", to.period="day"),
    1.2, seq(0, 25))

## qi = 500 bopd, Di = 3.91 nominal / year, b = 0.5, t = 5 years
hyperbolic.Np(rescale.by.time(500, from.period="day", to.period="year"),
    3.91, 0.5, 5)

## Di = 85% effective / year, b = 1.5, t = 6 months
hyperbolic.D(as.nominal(0.85), 1.5, 0.5)
}
