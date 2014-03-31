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

curtailed.q <- function (t.curtail, decline, t)
{
    q <- rep(decline$qi, length(t))
    q[t > t.curtail] <- arps.q(decline, t[t > t.curtail] - t.curtail)
    q
}

curtailed.Np <- function (t.curtail, decline, t)
{
    Np <- decline$qi * t
    Np[t > t.curtail] <- decline$qi * t.curtail +
      arps.Np(decline, t[t > t.curtail] - t.curtail)
    Np
}

curtailed.D <- function (t.curtail, decline, t)
{
    D <- rep(0, length(t))
    D[t > t.curtail] <- arps.D(decline, t[t > t.curtail] - t.curtail)
    D
}

curtail <- function (decline, t.curtail)
{
    res <- list(arps=decline, t.curtail=t.curtail)
    class(res) <- c("curtailed", "arps")
    res
}

arps.q.curtailed <- function(decl, t)
    curtailed.q(decl$t.curtail, decl$arps, t)

arps.Np.curtailed <- function(decl, t)
    curtailed.Np(decl$t.curtail, decl$arps, t)

arps.D.curtailed <- function(decl, t)
    curtailed.D(decl$t.curtail, decl$arps, t)

print.curtailed <- function(x, ...)
{
    cat("Curtailed ")
    print(x$arps)
    cat("with t_curtail = ", x$t.curtail, "\n", sep="")
    invisible(x)
}
