\name{format.arps}
\alias{format.arps}

\title{Format methods for Arps decline objects}

\description{
    Get human-readable representation of Arps decline-curve objects.
}

\usage{
\method{format}{arps}(x, \dots)
}

\arguments{
  \item{x}{Arps decline curve object as returned from \code{\link{arps.decline}}.}
  \item{\dots}{Arguments to additional \code{\link{format}} methods.}
}

\value{
    A character representation of \code{x}.
}

\seealso{
  \code{\link{format}}, \code{\link{print.arps}}, \code{\link{arps.decline}}.
}

\examples{
## exponential decline with
## qi = 1000 Mscf/d, Di = 95% effective / year
## rate for t from 0 to 25 days
decline <- arps.decline(1000,
    as.nominal(0.95, from.period="year", to.period="day"))
format(decline, digits=2)
}
