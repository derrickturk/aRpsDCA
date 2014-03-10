\name{aRpsDCA-package}
\alias{aRpsDCA-package}
\alias{aRpsDCA}

\docType{package}

\title{
Arps Decline Curve Analysis in R
}

\description{
Functions for Arps decline-curve analysis. Includes exponential, hyperbolic, harmonic, and hyperbolic-to-exponential models.
}

\details{
Index:
\tabular{ll}{
\code{\link{arps}}             \tab   Generic functions for Arps decline curves \cr
\code{\link{as.effective}}     \tab   Arps decline conversion from nominal to effective \cr
\code{\link{as.nominal}}       \tab   Arps decline conversion from effective to nominal \cr
\code{\link{bestfit}}          \tab   Best-fitting for Arps decline curves \cr
\code{\link{exponential}}      \tab   Arps exponential declines \cr
\code{\link{harmonic}}         \tab   Arps harmonic declines \cr
\code{\link{hyp2exp}}          \tab   Arps hyperbolic-to-exponential declines \cr
\code{\link{hyperbolic}}       \tab   Arps hyperbolic declines \cr
\code{\link{print.arps}}       \tab   Print representations of Arps decline curves \cr
\code{\link{rescale.by.time}}  \tab   Time unit conversion for DCA \cr
}
}

\author{
Derrick Turk | \bold{terminus data science, LLC} <\email{dwt@terminusdatascience.com}>
}

\references{
    SPEE REP#6 (\url{https://secure.spee.org/sites/default/files/wp-files/pdf/ReferencesResources/REP06-DeclineCurves.pdf})
}

\keyword{ package }

\examples{
    ## Plot semi-log rate-time and Cartesian rate-cumulative
    ## for a hyperbolic-to-exponential decline
    t <- seq(0, 10, 0.5)
    q <- hyp2exp.q(5000, as.nominal(0.90), 1.5, as.nominal(0.10), t)
    Np <- hyp2exp.Np(5000, as.nominal(0.90), 1.5, as.nominal(0.10), t)

    old.par <- par(ask=TRUE)
    plot(log(q) ~ t)
    plot(q ~ Np)
    par(old.par)
}
