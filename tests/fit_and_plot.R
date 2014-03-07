fitme.exponential.t <- seq(0, 5, 1 / 12) # 5 years
fitme.exponential.q <- exponential.q(
    1000, # Bbl/d
    as.nominal(0.70), # / year
    fitme.exponential.t
) * rnorm(n=length(fitme.exponential.t), mean=1, sd=0.1) # perturb

exponential.fit <- best.exponential(fitme.exponential.q, fitme.exponential.t)
cat(paste("SSE:", exponential.fit$sse))
plot(fitme.exponential.q ~ fitme.exponential.t, main="Exponential Fit",
     col="blue", log="y", xlab="Time", ylab="Rate")
lines(arps.q(exponential.fit$decline, fitme.exponential.t) ~ fitme.exponential.t,
      col="red")
legend("topright", pch=c(1, NA), lty=c(NA, 1), col=c("blue", "red"), legend=c("Actual", "Fit"))
