\name{curtailed}
\alias{curtailed}
\alias{curtail}
\alias{curtailed.q}
\alias{curtailed.Np}
\alias{curtailed.D}
\alias{arps.q.curtailed}
\alias{arps.Np.curtailed}
\alias{arps.D.curtailed}

\title{Arps decline curves with initial curtailment}

\description{
    Create decline curve objects and compute rates, cumulative production, and nominal declines for Arps decline curves with an initial period of curtailment to constant rate. The resulting decline curve models will impose a constant rate of production equal to the initial rate of the underlying model until the specified end-of-curtailment time.
}

\usage{
    curtail(decl, t.curtail)
    curtailed.q(decl, t.curtail, t)
    curtailed.Np(decl, t.curtail, t)
    curtailed.D(decl, t.curtail, t)
}

\arguments{
    \item{decl}{an Arps decline object as returned by \code{arps.decline}.}
    \item{t.curtail}{time to end of curtailment, in same units as \code{t} [time].}
    \item{t}{time at which to evaluate rate, cumulative, or nominal decline [time].}
}

\value{
    \code{curtail} returns an object having class \code{"arps"}, suitable
    for use as an argument to S3 methods listed in \code{\link{arps}}.

    \code{curtailed.q} returns the rate for each element of \code{t},
    under a decline described by \code{decl} and curtailed until time
    \code{t.curtail}.

    \code{curtailed.Np} returns the cumulative production for each element of
    \code{t}, under a decline described by \code{decl} and curtailed until time
    \code{t.curtail}.

    \code{curtailed.D} returns the nominal instantaneous decline for each
    element of \code{t}, under a decline described by \code{decl} and
    curtailed until time \code{t.curtail}.
}

\details{
    Assumes consistent units as described for underlying Arps decline types.
}

\seealso{
    \code{\link{arps}}, \code{\link{print.arps}}, \code{\link{exponential}},
    \code{\link{hyperbolic}}, \code{\link{hyp2exp}}.
}

\examples{
## qi = 1000 Mscf/d, Di = 95% effective / year, b = 1.2, t from 0 to 25 days,
## curtailed until 5 days
curtailed.q(arps.decline(
        1000,
        as.nominal(0.95, from.period="year", to.period="day"),
        1.2
    ),
    5, seq(0, 25))

## hyperbolic decline with
## qi = 500 bopd, Di = 3.91 nominal / year, b = 1.5,
## curtailed for 1 month
## cumulative production at t = 5 years
decline <- curtail(
    arps.decline(rescale.by.time(500, from="day", to="year"),
                 3.91, 1.5),
    (1 / 12)
)
arps.Np(decline, 5)
}
