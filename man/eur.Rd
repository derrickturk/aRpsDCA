\name{arps.eur}
\alias{arps.eur}
\alias{arps.t.el}

\title{EUR and time-to-limit for Arps decline curves}

\description{
    Evaluate estimated ultimate recovery and time to economic limit for Arps
    decline curve objects.
}

\usage{
arps.eur(decl, q.limit)
arps.t.el(decl, q.limit)
}

\arguments{
  \item{decl}{an Arps decline object as returned by \code{arps.decline} or \code{curtail}.}
  \item{q.limit}{economic limit rate [volume / time] in same units as \code{decl}.}
}

\value{
    \code{arps.eur} returns the total production for \code{decl} at the point in
    time when the economic limit \code{q.limit} is reached; that is,
    \code{arps.Np(decl, arps.t.el(decl, q.limit))}, in the same units as
    \code{q.limit}.

    \code{arps.t.el} returns the time until the economic limit \code{q.limit} is
    reached for decline \code{decl}.
}

\seealso{
  \code{\link{arps}}, \code{\link{curtail}}.
}

\examples{
## exponential decline with
## qi = 1000 Mscf/d, Di = 95% effective / year
## EUR with economic limit 10 Mscf/d
decline <- arps.decline(1000,
    as.nominal(0.95, from.period="year", to.period="day"))
arps.eur(decline, 10)

## hyperbolic decline with
## qi = 500 bopd, Di = 3.91 nominal / year, b = 1.5,
## time to reach economic limit of 3 bopd
decline <- arps.decline(rescale.by.time(500, from="day", to="year"),
    3.91, 1.5)
arps.t.el(decline, rescale.by.time(3, from="day", to="year"))
}

