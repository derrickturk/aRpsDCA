\name{rescale.by.time}
\alias{rescale.by.time}

\title{Time unit conversion for DCA}

\description{
    Scales rates, declines, and time periods from one time unit to another.
}

\usage{
    rescale.by.time(value,
        from.period=c("year", "month", "day"),
        to.period=c("year", "month", "day"),
        method=c("decline", "rate", "time"))
}

\arguments{
  \item{value}{rate [volume / time], Arps nominal decline [1 / time], or time to be rescaled.}
  \item{from.period}{time period for \code{value} (default \code{"year"}).}
  \item{to.period}{time period for result (default \code{"year"}).}
  \item{method}{scaling method to be applied, depending upon the type of \code{value} (default "decline").}
}

\value{
    Returns \code{value} scaled from \code{from.period} to \code{to.period} according to its type as specified by \code{method}.
}

\seealso{
  \code{\link{as.nominal}}, \code{\link{as.effective}}.
}

\examples{
    ## 3 MMscf/D to MMscf/year
    rescale.by.time(3, from.period="day", to.period="year", method="rate")

    ## Nominal decline of 3.2/year to nominal decline per month
    rescale.by.time(3.2, from.period="year", to.period="month", method="decline")

    ## 5 years in days
    rescale.by.time(5, from.period="year", to.period="month", method="time")
}
