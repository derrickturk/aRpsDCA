\name{arps.with.buildup}
\alias{arps.with.buildup}

\title{Arps declines with linear buildup period}

\description{
    Extend Arps decline curve objects by replacing early-time declines with a buildup period in which rate is a linear function of time.
}

\usage{
    arps.with.buildup(decl, initial.rate, time.to.peak)
}

\arguments{
    \item{decl}{an Arps decline object as produced by \code{arps.decline}.}
    \item{initial.rate}{initial rate [volume / time] (at time = 0) for buildup period.}
    \item{time.to.peak}{time to peak rate (i.e.~length of buildup period).}
}

\value{
    \code{arps.with.buildup} returns an object having class \code{"arps"},
    which may be used as an argument to methods such as \code{\link{arps.q}},
    \code{\link{arps.Np}}, \code{\link{arps.D}}, or \code{\link{print.arps}}.

    This object implements a decline curve which behaves as \code{decl} for all
    time greater than \code{time.to.peak}, but implements a linear buildup of
    rate interpolated between \code{initial.rate} at time zero and
    \code{arps.q(decl, time.to.peak)} at \code{time.to.peak}.
}

\seealso{
    \code{\link{arps.decline}}
}

\examples{
## hyperbolic decline with
## qi = 500 bopd, Di = 3.91 nominal / year, b = 1.5,
## cumulative production at t = 5 years
decline <- arps.decline(
    rescale.by.time(500, from="day", to="year", method="rate"),
    3.91, 1.5)
# add buildup from initial rate of 50 bopd, over 30 days
decline.with.buildup <- arps.with.buildup(decline,
    rescale.by.time(50, from="day", to="year", method="rate"),
    rescale.by.time(30, from="day", to="year", method="time"))
# forecast 5 years and compare
forecast.time <- seq(0, 5, 0.1)
plot(arps.q(decline, forecast.time) ~ forecast.time, log="y", type="l",
    lty="dashed", col="red")
lines(arps.q(decline.with.buildup, forecast.time) ~ forecast.time,
    lty="dotted", col="blue")
}
