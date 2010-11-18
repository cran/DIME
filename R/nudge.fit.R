nudge.fit <-
function(data, weights = NULL, pi = NULL, mu = NULL,
sigma = NULL, tol=1e-16, max.iter=2000, z = NULL)
{
  inudge.fit(data, K = 1, weights = weights, pi = pi, mu = mu,
    sigma = sigma, tol=tol, max.iter=max.iter, z = z);
}