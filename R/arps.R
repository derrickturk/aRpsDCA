HARMONIC_EPS <- 1e-3
EXPONENTIAL_EPS <- 1e-3

exponential.q <- function (qi, D, t)
{
    qi * exp(-D * t)
}

exponential.Np <- function (qi, D, t)
{
    qi / D * (1 - exp(-D * t))
}

harmonic.q <- function (qi, Di, t)
{
    qi / (1 + Di * t)
}

harmonic.Np <- function (qi, Di, t)
{
    qi / Di * log(1 + Di * t)
}

harmonic.D <- function (qi, Di, t)
{
    hyperbolic.D(qi, Di, 1, t)
}

hyperbolic.q <- function (qi, Di, b, t)
{
    if (abs(b - 1) < HARMONIC_EPS)
        harmonic.q(qi, Di, t)
    else if (abs(b) < EXPONENTIAL_EPS)
        exponential.q(qi, Di, t)
    else
        qi * (1 + b * Di * t) ^ (-1/b)
}

hyperbolic.Np <- function (qi, Di, b, t)
{
    if (abs(b - 1) < HARMONIC_EPS)
        harmonic.Np(qi, Di, t)
    else if (abs(b) < EXPONENTIAL_EPS)
        exponential.Np(qi, Di, t)
    else
        (qi / ((1 - b) * Di)) * (1 - (1 + b * Di * t) ^ (1 - (1/b)))
}

hyperbolic.D <- function (qi, Di, b, t)
{
    Di / (1 + b * Di * t)
}

hyp2exp.transition <- function (qi, Di, b, Df)
{
    (Di / Df - 1) / (b * Di)
}

hyp2exp.q <- function (qi, Di, b, Df, t)
{
    t.trans <- hyp2exp.transition(qi, Di, b, Df)
    q.trans <- hyperbolic.q(qi, Di, b, t.trans)

    q <- hyperbolic.q(qi, Di, b, t)
    q[t > t.trans] <- exponential.q(q.trans, Df, t[t > t.trans] - t.trans)

    q
}

hyp2exp.Np <- function (qi, Di, b, Df, t)
{
    t.trans <- hyp2exp.transition(qi, Di, b, Df)
    q.trans <- hyperbolic.q(qi, Di, b, t.trans)
    Np.trans <- hyperbolic.Np(qi, Di, b, t.trans)

    Np <- hyperbolic.Np(qi, Di, b, t)
    Np[t > t.trans] <- Np.trans
        + exponential.Np(q.trans, Df, t[t > t.trans] - t.trans)

    Np
}

hyp2exp.D <- function (qi, Di, b, Df, t)
{
    t.trans <- hyp2exp.transition(qi, Di, b, Df)
    D <- hyperbolic.D(qi, Di, b, t)
    D[t > t.trans] <- Df

    D
}
