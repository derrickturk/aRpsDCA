# aRpsDCA
# Copyright (C) 2014 dwt | terminus data science, LLC
# <dwt [at] terminusdatascience.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA

arps.decline <- function(qi, Di, b=NA, Df=NA)
{
    if (is.na(b) && !is.na(Df))
        stop("Must specify b for hyp2exp decline.")

    if (is.na(b)) {
        res <- list(qi=qi, D=Di)
        class(res) <- c("exponential", "arps")
    } else if (is.na(Df)) {
        res <- list(qi=qi, Di=Di, b=b)
        class(res) <- c("hyperbolic", "arps")
    } else {
        res <- list(qi=qi, Di=Di, b=b, Df=Df)
        class(res) <- c("hyp2exp", "arps")
    }

    res
}

arps.q <- function(decl, t)
{
    UseMethod("arps.q")
}

arps.Np <- function(decl, t)
{
    UseMethod("arps.Np")
}

arps.D <- function(decl, t)
{
    UseMethod("arps.D")
}

arps.q.arps <- function(decl, t) stop("Need specific decline class.")
arps.Np.arps <- function(decl, t) stop("Need specific decline class.")
arps.D.arps <- function(decl, t) stop("Need specific decline class.")

arps.q.exponential <- function(decl, t) exponential.q(decl$qi, decl$D, t)
arps.q.hyperbolic <- function(decl, t) hyperbolic.q(decl$qi, decl$Di, decl$b, t)
arps.q.hyp2exp <- function(decl, t) hyp2exp.q(decl$qi, decl$Di, decl$b, decl$Df, t)

arps.Np.exponential <- function(decl, t) exponential.Np(decl$qi, decl$D, t)
arps.Np.hyperbolic <- function(decl, t) hyperbolic.Np(decl$qi, decl$Di, decl$b, t)
arps.Np.hyp2exp <- function(decl, t) hyp2exp.Np(decl$qi, decl$Di, decl$b, decl$Df, t)

arps.D.exponential <- function(decl, t) decl$D
arps.D.hyperbolic <- function(decl, t) hyperbolic.D(decl$Di, decl$b, t)
arps.D.hyp2exp <- function(decl, t) hyp2exp.D(decl$Di, decl$b, decl$Df, t)

format.arps <- function(x, ...)
{
    paste("Arps decline:", format(unclass(x), ...), sep="\n")
}

format.exponential <- function(x, ...)
{
    paste("Arps exponential decline: <qi = ",
          format(x$qi, ...),
          ", D = ",
          format(x$D, ...),
          ">",
          sep="")
}

format.hyperbolic <- function(x, ...)
{
    paste("Arps exponential decline: <qi = ",
          format(x$qi, ...),
          ", Di = ",
          format(x$Di, ...),
          ", b = ",
          format(x$b, ...),
          ">",
          sep="")
}

format.hyp2exp <- function(x, ...)
{
    paste("Arps exponential decline: <qi = ",
          format(x$qi, ...),
          ", Di = ",
          format(x$Di, ...),
          ", b = ",
          format(x$b, ...),
          ", Df = ",
          format(x$Df, ...),
          ">",
          sep="")
}

print.arps <- function(x, ...)
{
    print(format(x, ...))
    invisible(x)
}
