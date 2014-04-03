\name{print.arps}
\alias{print.arps}

\title{Print methods for Arps decline objects}

\description{
    Print human-readable representation of Arps decline-curve objects using \code{\link{format.arps}}.
}

\usage{
\method{print}{arps}(x, \dots)
}

\arguments{
  \item{x}{Arps decline curve object as returned from \code{\link{arps.decline}}.}
  \item{\dots}{Arguments to \code{\link{format.arps}}.}
}

\value{
    Invisibly (see \code{\link{invisible}}) returns \code{x}.
}

\seealso{
  \code{\link{print}}, \code{\link{format.arps}}, \code{\link{arps.decline}}.
}

\examples{
## exponential decline with
## qi = 1000 Mscf/d, Di = 95% effective / year
## rate for t from 0 to 25 days
decline <- arps.decline(1000,
    as.nominal(0.95, from.period="year", to.period="day"))
print(decline)
}
