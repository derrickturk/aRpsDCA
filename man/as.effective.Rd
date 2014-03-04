\name{as.effective}
\alias{as.effective}

\title{Arps decline conversion from nominal to effective}

\description{
    Convert nominal to effective Arps decline.
}

\usage{
    as.effective(D.nom,
        from.period=c("year", "month", "day"),
        to.period=c("year", "month", "day"))
}

\arguments{
  \item{D.nom}{nominal Arps decline [1 / time].}
  \item{from.period}{time period for \code{D.nom} (default \code{"year"}).}
  \item{to.period}{time period for result (default \code{"year"}).}
}

\value{
    Returns the tangent effective Arps decline rate in units of [1 / \code{to.period}].
}

\details{
    The result should be applied as a tangent effective decline (see SPEE REP#6 [\link{https://secure.spee.org/sites/default/files/wp-files/pdf/ReferencesResources/REP06-DeclineCurves.pdf}]) specified in fractional (i.e. 95\% = 0.95) units.

    When appropriate, internally uses \code{\link{rescale.by.time}} to perform time unit conversion.
}

\seealso{
  \code{\link{as.nominal}}, \code{\link{rescale.by.time}}.
}

\examples{
    ## 0.008 nominal daily decline to tangent effective yearly decline
    as.effective(0.008, from.period="day", to.period="year")
}
