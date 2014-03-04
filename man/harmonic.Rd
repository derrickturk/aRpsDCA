\name{harmonic}
\alias{harmonic}
\alias{harmonic.q}
\alias{harmonic.Np}
\alias{harmonic.D}

\title{Arps harmonic declines}

\description{
    Compute rates, cumulative production values, and instantaneous nominal declines for Arps harmonic decline curves (i.e. hyperbolic with b = 1).
}

\usage{
    harmonic.q(qi, Di, t)
    harmonic.Np(qi, Di, t)
    harmonic.D(Di, t)
}

\arguments{
  \item{qi}{initial rate [volume / time], i.e. q(t = 0).}
  \item{Di}{initial nominal Arps decline exponent [1 / time].}
  \item{t}{time at which to evaluate rate or cumulative [time].}
}

\value{
    \code{harmonic.q} returns the rate for each element of \code{t},
    in the same units as \code{qi}.

    \code{harmonic.Np} returns the cumulative production for each element of
    \code{t}, in the same units as \code{qi * t}.

    \code{harmonic.D} returns the nominal instantaneous decline for each
    element of \code{t}. This can be converted to effective decline and
    rescaled in time by use of \code{\link{as.effective}} and
    \code{\link{rescale.by.time}}.
}

\details{
    Assumes consistent units of time between \code{qi}, \code{D}, and \code{t}. To convert, see the decline-rate conversion functions referenced below.
}

\seealso{
  \code{\link{as.effective}}, \code{\link{as.nominal}}, \code{\link{rescale.by.time}}.
}

\examples{
    ## qi = 1000 Mscf/d, Di = 95% effective / year, t from 0 to 25 days
    harmonic.q(1000, as.nominal(0.95, from.period="year", to.period="day"), seq(0, 25))

    ## qi = 500 bopd, Di = 3.91 nominal / year, t = 5 years
    harmonic.Np(rescale.by.time(500, from.period="day", to.period="year"), 3.91, 5)

    ## Di = 85% effective / year, t = 6 months
    harmonic.D(as.nominal(0.85), 0.5)
}

