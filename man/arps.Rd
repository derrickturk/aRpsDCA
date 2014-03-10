\name{arps}
\alias{arps}
\alias{arps.decline}
\alias{arps.q}
\alias{arps.Np}
\alias{arps.D}

\title{Arps decline classes and S3 methods}

\description{
    Create Arps decline curve objects and compute rates, cumulative production, and nominal declines.
}

\usage{
    arps.decline(qi, Di, b=NA, Df=NA)
    arps.q(decl, t)
    arps.Np(decl, t)
    arps.D(decl, t)
}

\arguments{
  \item{qi}{initial rate [volume / time], i.e. q(t = 0).}
  \item{Di}{nominal Arps decline exponent [1 / time].}
  \item{b}{Arps hyperbolic exponent.}
  \item{Df}{nominal Arps decline exponent [1 / time].}
  \item{t}{time at which to evaluate rate, cumulative, or nominal decline [time].}
  \item{decl}{An Arps decline object as returned by \code{arps.decline}.}
}

\value{
    \code{arps.decline} returns an object having class \code{"arps"}, suitable
    for use as an argument to S3 methods discussed here.

    \code{q.arps} returns the rate for each element of \code{t} applying
    decline \code{decl}, in the same units as the value of \code{qi} for
    \code{decl}.

    \code{Np.arps} returns the cumulative production for each element of \code{t}
    applying decline \code{decl}, in the same units as the value of
    \code{qi * t} for \code{decl}.

    \code{D.arps} returns the nominal decline for each element of \code{t}
    applying decline \code{decl}, in the same units as the value of \code{Di}
    for \code{decl}.
}

\details{
    Depending on whether arguments \code{b} and \code{Df} are supplied,
    \code{arps.decline} will select an exponential, hyperbolic, or
    hyperbolic-to-exponential decline and return an object appropriately.
    The returned object will have class \code{"exponential"},
    \code{"hyperbolic"}, or \code{"hyp2exp"} in addition to class
    \code{"arps"}.

    Assumes consistent units of time between \code{qi}, \code{Di}, \code{Df},
    and \code{t}.
    To convert, see the decline-rate conversion functions referenced below.
}

\seealso{
  \code{\link{print.arps}}, \code{\link{exponential}}, \code{\link{hyperbolic}},
  \code{\link{hyp2exp}}, \code{\link{as.effective}}, \code{\link{as.nominal}},
  \code{\link{rescale.by.time}}.
}

\examples{
    ## exponential decline with
    ## qi = 1000 Mscf/d, Di = 95% effective / year
    ## rate for t from 0 to 25 days
    decline <- arps.decline(1000,
        as.nominal(0.95, from.period="year", to.period="day"))
    arps.q(decline, seq(0, 25))

    ## hyperbolic decline with
    ## qi = 500 bopd, Di = 3.91 nominal / year, b = 1.5,
    ## cumulative production at t = 5 years
    decline <- arps.decline(rescale.by.time(500, from="day", to="year"),
        3.91, 1.5)
    arps.Np(decline, 5)
}
