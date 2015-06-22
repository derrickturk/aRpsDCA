# aRpsDCA
# Copyright (C) 2015 dwt | terminus data science, LLC
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
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1])),
                         # Di = decline from first to second point

                    # cost function
                 function (guess) sse(q, exponential.q(guess[1], guess[2], t)),

                 lower=c( # lower bounds
                   0, # qi > 0
                   0), # D > 0

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10) # = 0.99995 / [time] effective
    )

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2]),
         sse=res$objective)
}

best.hyperbolic <- function(q, t)
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                    # cost function
                 function (guess)
                     sse(q, hyperbolic.q(guess[1], guess[2], guess[3], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0,  # Di > 0
                   0), # b > 0

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   5)  # don't get crazy
    )

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
         sse=res$objective)
}

best.hyp2exp <- function(q, t)
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                    # cost function
                 function (guess)
                     sse(q,
                         hyp2exp.q(guess[1], guess[2], guess[3], guess[4], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0.35,  # Di > 0
                   0,  # b > 0
                   0), # Df > 0

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   5,  # don't get crazy
                   0.35) # likewise
    )

    list(decline=arps.decline(qi=res$par[1],
                              Di=res$par[2],
                              b=res$par[3],
                              Df=res$par[4]),
         sse=res$objective)
}

best.exponential.curtailed <- function(q, t)
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   t[2]  # t.curtail = second t in vector
                   ),

                    # cost function
                 function (guess)
                     sse(q,
                         curtailed.q(arps.decline(guess[1], guess[2]),
                                     guess[3], t)),

                 lower=c( # lower bounds
                   0, # qi > 0
                   0, # D > 0
                   0  # t.curtail > 0
                 ),

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   t[length(t)]
                 )
    )

    list(decline=curtail(arps.decline(qi=res$par[1], Di=res$par[2]),
                         res$par[3]),
         sse=res$objective)
}

best.hyperbolic.curtailed <- function(q, t)
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5,  # right-ish for a lot of wells currently coming on
                   t[2]  # t.curtail = second t in vector
                   ),

                    # cost function
                 function (guess)
                     sse(q,
                         curtailed.q(
                                 arps.decline(guess[1], guess[2], guess[3]),
                                 guess[4], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0,  # Di > 0
                   0,  # b > 0
                   0   # t.curtail > 0
                 ),

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   5,  # don't get crazy
                   t[length(t)]
                 )
    )

    list(decline=curtail(
             arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
             res$par[4]),
         sse=res$objective)
}

best.hyp2exp.curtailed <- function(q, t)
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   (log(q[2]) - log(q[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1,  # Df = about 9% effective
                   t[2]
                   ),

                    # cost function
                 function (guess)
                     sse(q,
                         curtailed.q(
                             arps.decline(guess[1], guess[2], guess[3], guess[4]),
                             guess[5], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0.35,  # Di > 0
                   0,  # b > 0
                   0,  # Df > 0
                   0
                 ),

                 upper=c( # upper bounds
                   max(q) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   5,  # don't get crazy
                   0.35, # likewise
                   t[length(t)]
                 )
    )

    list(decline=curtail(arps.decline(qi=res$par[1],
                                      Di=res$par[2],
                                      b=res$par[3],
                                      Df=res$par[4]),
                         res$par[5]),
         sse=res$objective)
}

best.fit <- function(q, t)
{
    exp <- best.exponential(q, t)
    hyp <- best.hyperbolic(q, t)
    h2e <- best.hyp2exp(q, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.curtailed.fit <- function(q, t)
{
    exp <- best.exponential.curtailed(q, t)
    hyp <- best.hyperbolic.curtailed(q, t)
    h2e <- best.hyp2exp.curtailed(q, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.exponential.from.Np <- function(Np, t)
{
    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1], # qi = Np(t = first t in vector)
                   (log(Np[2] - Np[1]) - log(Np[1])) / (t[2] - t[1])),
                         # Di = decline from first to second point

                    # cost function
                 function (guess) sse(Np, exponential.Np(guess[1], guess[2], t)),

                 lower=c( # lower bounds
                   0, # qi > 0
                   0), # D > 0

                 upper=c( # upper bounds
                   max(Np) * 3, # qi < max(Np) * 3
                   10) # = 0.99995 / [time] effective
    )

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2]),
         sse=res$objective)
}

best.hyperbolic.from.Np <- function(Np, t)
{
    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1], # qi = q(t = first t in vector)
                   (log(Np[2] - Np[1]) - log(Np[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                    # cost function
                 function (guess)
                     sse(Np, hyperbolic.Np(guess[1], guess[2], guess[3], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0,  # Di > 0
                   0), # b > 0

                 upper=c( # upper bounds
                   max(Np) * 3, # qi < max(Np) * 3
                   10, # = 0.99995 / [time] effective
                   5)  # don't get crazy
    )

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
         sse=res$objective)
}

best.hyp2exp.from.Np <- function(Np, t)
{
    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1], # qi = Np(t = first t in vector)
                   (log(Np[2] - Np[1]) - log(Np[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                    # cost function
                 function (guess)
                     sse(Np,
                         hyp2exp.Np(guess[1], guess[2], guess[3], guess[4], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0.35,  # Di > 0
                   0,  # b > 0
                   0), # Df > 0

                 upper=c( # upper bounds
                   max(Np) * 3, # qi < max(Np) * 3
                   10, # = 0.99995 / [time] effective
                   5,  # don't get crazy
                   0.35) # likewise
    )

    list(decline=arps.decline(qi=res$par[1],
                              Di=res$par[2],
                              b=res$par[3],
                              Df=res$par[4]),
         sse=res$objective)
}

best.exponential.curtailed.from.Np <- function(Np, t)
{
    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1], # qi = Np(t = first t in vector)
                   (log(Np[2]) - log(Np[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   t[2]  # t.curtail = second t in vector
                   ),

                    # cost function
                 function (guess)
                     sse(Np,
                         curtailed.Np(arps.decline(guess[1], guess[2]),
                                     guess[3], t)),

                 lower=c( # lower bounds
                   0, # qi > 0
                   0, # D > 0
                   0  # t.curtail > 0
                 ),

                 upper=c( # upper bounds
                   max(Np) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   t[length(t)]
                 )
    )

    list(decline=curtail(arps.decline(qi=res$par[1], Di=res$par[2]),
                         res$par[3]),
         sse=res$objective)
}

best.hyperbolic.curtailed.from.Np <- function(Np, t)
{
    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1], # Np = Np(t = first t in vector)
                   (log(Np[2]) - log(Np[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5,  # right-ish for a lot of wells currently coming on
                   t[2]  # t.curtail = second t in vector
                   ),

                    # cost function
                 function (guess)
                     sse(Np,
                         curtailed.Np(
                                 arps.decline(guess[1], guess[2], guess[3]),
                                 guess[4], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0,  # Di > 0
                   0,  # b > 0
                   0   # t.curtail > 0
                 ),

                 upper=c( # upper bounds
                   max(Np) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   5,  # don't get crazy
                   t[length(t)]
                 )
    )

    list(decline=curtail(
             arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
             res$par[4]),
         sse=res$objective)
}

best.hyp2exp.curtailed.from.Np <- function(Np, t)
{
    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1], # Np = Np(t = first t in vector)
                   (log(Np[2]) - log(Np[1])) / (t[2] - t[1]),
                         # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1,  # Df = about 9% effective
                   t[2]
                   ),

                    # cost function
                 function (guess)
                     sse(Np,
                         curtailed.Np(
                             arps.decline(guess[1], guess[2], guess[3], guess[4]),
                             guess[5], t)),

                 lower=c( # lower bounds
                   0,  # qi > 0
                   0.35,  # Di > 0
                   0,  # b > 0
                   0,  # Df > 0
                   0
                 ),

                 upper=c( # upper bounds
                   max(Np) * 3, # qi < qmax * 3
                   10, # = 0.99995 / [time] effective
                   5,  # don't get crazy
                   0.35, # likewise
                   t[length(t)]
                 )
    )

    list(decline=curtail(arps.decline(qi=res$par[1],
                                      Di=res$par[2],
                                      b=res$par[3],
                                      Df=res$par[4]),
                         res$par[5]),
         sse=res$objective)
}

best.fit.from.Np <- function(Np, t)
{
    exp <- best.exponential.from.Np(Np, t)
    hyp <- best.hyperbolic.from.Np(Np, t)
    h2e <- best.hyp2exp.from.Np(Np, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.curtailed.fit.from.Np <- function(Np, t)
{
    exp <- best.exponential.curtailed.from.Np(Np, t)
    hyp <- best.hyperbolic.curtailed.from.Np(Np, t)
    h2e <- best.hyp2exp.curtailed.from.Np(Np, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}
