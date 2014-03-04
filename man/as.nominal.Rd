\name{as.nominal}
\alias{as.nominal}

\title{Arps decline conversion from effective to nominal}

\description{
    Convert effective to nominal Arps decline.
}

\usage{
    as.nominal(D.eff,
        from.period=c("year", "month", "day"),
        to.period=c("year", "month", "day"))
}

\arguments{
  \item{D.eff}{tangent effective Arps decline [1 / time].}
  \item{from.period}{time period for \code{D.eff} (default \code{"year"}).}
  \item{to.period}{time period for result (default \code{"year"}).}
}

\value{
    Returns the Arps nominal decline rate in units of [1 / \code{to.period}].
}

\details{
    \code{D.eff} should be specified in fractional (i.e. 95\% = 0.95) units as a tangent effective decline (see SPEE REP#6 [\link{https://secure.spee.org/sites/default/files/wp-files/pdf/ReferencesResources/REP06-DeclineCurves.pdf}]).

    When appropriate, internally uses \code{\link{rescale.by.time}} to perform time unit conversion.
}

\seealso{
  \code{\link{as.effective}}, \code{\link{rescale.by.time}}.
}

\examples{
    ## 95% / year effective decline to nominal daily decline
    as.nominal(0.95, from.period="year", to.period="day")
}
