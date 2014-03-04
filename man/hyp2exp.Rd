\name{hyp2exp}
\alias{hyp2exp}
\alias{hyp2exp.q}
\alias{hyp2exp.Np}
\alias{hyp2exp.D}
\alias{hyp2exp.transition}

\title{Arps hyperbolic-to-exponential declines}

\description{
    Compute rates, cumulative production values, instantaneous nominal declines, and transition times for Arps hyperbolic-to-exponential decline curves.
}

\usage{
    hyp2exp.q(qi, Di, b, Df, t)
    hyp2exp.Np(qi, Di, b, Df, t)
    hyp2exp.D(qi, Di, b, Df, t)
    hyp2exp.transition(qi, Di, b, Df)
}

\arguments{
  \item{qi}{initial rate [volume / time], i.e. q(t = 0).}
  \item{Di}{initial nominal Arps decline exponent [1 / time].}
  \item{b}{Arps hyperbolic exponent.}
  \item{Df}{final nominal Arps decline exponent [1 / time].}
  \item{t}{time at which to evaluate rate or cumulative [time].}
}

\value{
    \code{hyp2exp.q} returns the rate for each element of \code{t},
    in the same units as \code{qi}.

    \code{hyp2exp.Np} returns the cumulative production for each element of
    \code{t}, in the same units as \code{qi * t}.

    \code{hyp2exp.D} returns the nominal instantaneous decline for each
    element of \code{t}. This can be converted to effective decline and
    rescaled in time by use of \code{\link{as.effective}} and
    \code{\link{rescale.by.time}}.

    \code{hyp2exp.transition} returns the transition time (from hyperbolic to
    exponential decline), in the same units as \code{1 / Di}.
}

\details{
    Assumes consistent units of time between \code{qi}, \code{Di}, \code{Df}, and \code{t}. To convert, see the decline-rate conversion functions referenced below.

    When appropriate, internally uses \code{\link{harmonic.q}} and \code{\link{harmonic.Np}} to avoid singularities in calculations for \code{b} near 1.
}

\seealso{
  \code{\link{as.effective}}, \code{\link{as.nominal}}, \code{\link{rescale.by.time}}.
}

\examples{
    ## qi = 1000 Mscf/d, Di = 95% effective / year, b = 1.2,
    ## Df = 15% effective / year, t from 0 to 25 years
    ## result in Mscf/d
    hyp2exp.q(1000, as.nominal(0.95), 1.2, as.nominal(0.15), seq(0, 25))

    ## qi = 500 bopd, Di = 3.91 nominal / year, b = 0.5,
    ## Df = 0.3 nominal / year, t = 20 years
    hyp2exp.Np(rescale.by.time(500, from.period="day", to.period="year"),
        3.91, 0.5, 0.3, 20)

    ## Di = 85% effective / year, b = 1.5, Df = 15% effective / year,
    ## t = 6 and 48 months
    hyp2exp.D(as.nominal(0.85), 1.5, c(0.5, 4))

    ## Di = 85% effective / year, b = 1.5, Df = 5% effective / year
    hyp2exp.transition(as.nominal(0.85), 1.5, as.nominal(0.05))
}
