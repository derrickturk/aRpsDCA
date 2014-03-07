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

arps.q.exponential <- function(decl, t) do.call(exponential.q, c(decl, list(t=t)))
arps.q.hyperbolic <- function(decl, t) do.call(hyperbolic.q, c(decl, list(t=t)))
arps.q.hyp2exp <- function(decl, t) do.call(hyp2exp.q, c(decl, list(t=t)))

arps.Np.exponential <- function(decl, t) do.call(exponential.Np, c(decl, list(t=t)))
arps.Np.hyperbolic <- function(decl, t) do.call(hyperbolic.Np, c(decl, list(t=t)))
arps.Np.hyp2exp <- function(decl, t) do.call(hyp2exp.Np, c(decl, list(t=t)))

arps.D.exponential <- function(decl, t) decl$D
arps.D.hyperbolic <- function(decl, t) do.call(hyperbolic.D, c(decl, list(t=t)))
arps.D.hyp2exp <- function(decl, t) do.call(hyp2exp.D, c(decl, list(t=t)))
