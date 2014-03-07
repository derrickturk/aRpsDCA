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

sse <- function(q, forecast)
{
    sum((q - forecast) ^ 2)
}

best.exponential <- function(q, t)
{
    if (length(q) != length(t) || length(q) <= 1)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1])),
                         # Di = decline from first to second point

                    # cost function
                 function (guess) sse(q, exponential.q(guess[1], guess[2], t)),

                 lower=c( # lower bounds
                   0, # qi > 0
                   0), # qi < qmax * 3

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10) # = 0.99995 / [time] effective
    )

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2]),
         sse=res$objective)
}
