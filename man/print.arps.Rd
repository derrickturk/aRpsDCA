\name{print.arps}
\alias{print}
\alias{print.arps}

\title{Print methods for Arps decline objects}

\description{
    Print human-readable representation of Arps decline-curve objects.
}

\usage{
    print(x, ...)
}

\arguments{
  \item{x}{Arps decline curve object as returned from \code{\link{arps.decline}}.}
  \item{...}{Arguments to additional methods.}
}

\value{
    Invisibly (see \code{\link{invisible}}) returns \code{x}.
}

\seealso{
  \code{\link{base::print}}, \code{\link{arps.decline}}.
}

\examples{
    ## exponential decline with
    ## qi = 1000 Mscf/d, Di = 95% effective / year
    ## rate for t from 0 to 25 days
    decline <- arps.decline(1000,
        as.nominal(0.95, from.period="year", to.period="day"))
    print(decline)
}
