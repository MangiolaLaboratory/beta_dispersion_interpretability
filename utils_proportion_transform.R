## Utilities: proportion + variance-stabilizing transforms for compositional counts
##
## These helpers are intentionally dependency-free (base R only) so they can be
## sourced from scripts and Quarto documents.
##
## Conventions used in this repo:
## - X_ij (count) with N_i (library_size)
## - p_hat = X_ij / N_i
## - arcsin-sqrt transform: asin(sqrt(p))

.clamp01 <- function(x, eps = 0) {
  x <- pmax(x, 0 + eps)
  x <- pmin(x, 1 - eps)
  x
}

#' Convert counts to proportions p_hat = count / library_size
#'
#' @param count numeric vector (>= 0)
#' @param library_size numeric vector (> 0), recycled to length(count) if needed
#' @param eps numeric scalar. If > 0, proportions are clamped to [eps, 1-eps].
#'
#' @return numeric vector of proportions in [0, 1] (or [eps, 1-eps] if eps > 0)
count_to_proportion <- function(count, library_size, eps = 0) {
  if (any(library_size <= 0, na.rm = TRUE)) {
    stop("library_size must be > 0 for all entries.")
  }
  p_hat <- count / library_size
  if (eps > 0) p_hat <- .clamp01(p_hat, eps = eps)
  p_hat
}

#' Arcsin-sqrt variance-stabilizing transform for proportions
#'
#' @param p numeric vector of proportions in [0, 1]
#' @param eps numeric scalar. If > 0, p is clamped to [eps, 1-eps] before transform.
#'
#' @return numeric vector asin(sqrt(p))
arcsin_sqrt <- function(p, eps = 0) {
  if (eps > 0) p <- .clamp01(p, eps = eps)
  asin(sqrt(p))
}

#' Residual (generic): subtraction on any scale
#'
#' This is intentionally transformation-agnostic. Whatever transformation you want
#' (e.g., proportions, arcsin-sqrt, CLR, etc.) should be applied BEFORE calling this.
#'
#' @param observed numeric vector (already on the chosen scale)
#' @param expected numeric vector (already on the same scale as observed)
#'
#' @return observed - expected
residual_subtract <- function(observed, expected) {
  observed - expected
}

#' Adjust residuals for a mean–dispersion association in unconstrained space
#'
#' We assume a (slope-only) association in unconstrained log-linear predictor space:
#' \deqn{\log(\sigma_{\text{rel}}) = -k \cdot \mu}
#'
#' and scale residuals by:
#' \deqn{r^* = r / \sqrt{\sigma_{\text{rel}}}}
#'
#' where:
#' - \eqn{r} is a residual on any chosen scale (e.g. arcsin–sqrt),
#' - \eqn{\mu} is the unconstrained log-linear predictor ("mu_inv_softmax"),
#' - \eqn{k} is `log_dispersion_assoc`.
#'
#' This is intentionally **independent** of how residuals were constructed.
#'
#' @param residual numeric vector/matrix of residuals (any scale).
#' @param mu_inv_softmax numeric vector/matrix of unconstrained log-linear predictors.
#'   Must be conformable with `residual` (same length or same dimensions).
#' @param log_dispersion_assoc numeric scalar `k` in log(sigma) = -k * mu + c.
#'   Only the slope is used here (intercept intentionally excluded).
#' @param clamp_log_sigma_rel optional numeric length-2 vector giving bounds for
#'   log(sigma_rel) to prevent overflow (e.g. c(-30, 30)). Default NULL = no clamp.
#'
#' @return residual adjusted for the mean–dispersion association.
adjust_residual_assoc <- function(
  residual,
  mu_inv_softmax,
  log_dispersion_assoc,
  clamp_log_sigma_rel = NULL
) {
  if (!is.numeric(residual)) stop("residual must be numeric.")
  if (!is.numeric(mu_inv_softmax)) stop("mu_inv_softmax must be numeric.")
  if (!is.numeric(log_dispersion_assoc) || length(log_dispersion_assoc) != 1) {
    stop("log_dispersion_assoc must be a numeric scalar.")
  }

  log_sigma_rel <- -log_dispersion_assoc * mu_inv_softmax

  if (!is.null(clamp_log_sigma_rel)) {
    if (!is.numeric(clamp_log_sigma_rel) || length(clamp_log_sigma_rel) != 2) {
      stop("clamp_log_sigma_rel must be a numeric vector of length 2, e.g. c(-30, 30).")
    }
    log_sigma_rel <- pmax(log_sigma_rel, clamp_log_sigma_rel[1])
    log_sigma_rel <- pmin(log_sigma_rel, clamp_log_sigma_rel[2])
  }

  sigma_rel <- exp(log_sigma_rel)
  residual / sqrt(sigma_rel)
}

#' Translate residuals back to a group/taxon expected location
#'
#' This is the "add the centroid back" step:
#' \deqn{y = r + c}
#'
#' where:
#' - \eqn{r} is a residual on some scale (e.g. arcsin–sqrt),
#' - \eqn{c} is the expected abundance/location on the *same scale*.
#'
#' Note: the expected abundances should be transformed externally with the same
#' transform used for the observed proportions (e.g. `arcsin_sqrt()`), so this
#' function remains transformation-agnostic.
#'
#' @param residual numeric vector/matrix of residuals.
#' @param expected_location numeric vector/matrix of expected locations on the same scale.
#'   Must be conformable with `residual` (same length or same dimensions).
#'
#' @return translated values `residual + expected_location`.
translate_residual_location <- function(residual, expected_location) {
  if (!is.numeric(residual)) stop("residual must be numeric.")
  if (!is.numeric(expected_location)) stop("expected_location must be numeric.")
  residual + expected_location
}

#' Convert dispersion (sigma) to beta-binomial ICC parameter rho
#'
#' In this repo, `sigma` denotes dispersion (higher = more overdispersion).
#' A convenient reparameterization is:
#' \deqn{\rho = \sigma / (1 + \sigma)}
#'
#' @param sigma numeric vector/matrix of dispersion values (> 0)
#' @return rho in (0, 1)
sigma_to_rho <- function(sigma) {
  if (!is.numeric(sigma)) stop("sigma must be numeric.")
  if (any(sigma <= 0, na.rm = TRUE)) stop("sigma must be > 0.")
  # Numerically stable for very large sigma (including Inf):
  # sigma/(1+sigma) can become Inf/Inf -> NaN when sigma is Inf.
  rho <- sigma / (1 + sigma)
  rho[is.infinite(sigma)] <- 1

  # Clamp away from exact boundaries to keep downstream variance formulas well-defined.
  eps <- 1e-12
  rho <- pmax(pmin(rho, 1 - eps), eps)
  rho
}

#' Approximate SD on arcsin-sqrt scale under beta-binomial sampling
#'
#' For proportions \eqn{\hat p = X/N} with beta-binomial overdispersion \eqn{\rho},
#' a delta-method approximation gives:
#' \deqn{\mathrm{Var}(\arcsin\sqrt{\hat p}) \approx (1 + (N-1)\rho)/(4N)}
#'
#' @param library_size numeric vector of N (> 0)
#' @param rho numeric vector of rho in (0, 1)
#' @return numeric vector SD on arcsin-sqrt scale
bb_arcsin_sqrt_sd <- function(library_size, rho) {
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(rho <= 0 | rho >= 1, na.rm = TRUE)) stop("rho must be in (0, 1).")

  sqrt((1 + (library_size - 1) * rho) / (4 * library_size))
}

#' Beta-binomial variance of the sample proportion p_hat = X/N
#'
#' Using rho parameterization:
#' \deqn{Var(p_hat) = pi(1-pi) * (1 + (N-1)rho) / N}
#'
#' @param pi numeric vector in [0,1]
#' @param library_size numeric vector of N (>0)
#' @param rho numeric vector in (0,1)
#' @param eps clamp for pi away from 0/1
#'
#' @return numeric vector variance of p_hat
bb_var_p_hat <- function(pi, library_size, rho, eps = 1e-12) {
  if (!is.numeric(pi)) stop("pi must be numeric.")
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")
  if (any(rho <= 0 | rho >= 1, na.rm = TRUE)) stop("rho must be in (0, 1).")

  pi <- .clamp01(pi, eps = eps)
  pi * (1 - pi) * (1 + (library_size - 1) * rho) / library_size
}

#' Posterior mean of the underlying proportion p under a beta-binomial model
#'
#' Model:
#' - p ~ Beta(alpha, beta)
#' - X | p ~ Binomial(N, p)
#'
#' Using the common beta-binomial ICC parameterization:
#' \deqn{\rho = 1 / (\alpha + \beta + 1)}
#'
#' Given mean \eqn{\pi} and \eqn{\rho} we have:
#' \deqn{\phi = \alpha + \beta = 1/\rho - 1,\ \alpha = \pi\phi,\ \beta = (1-\pi)\phi}
#'
#' Posterior:
#' \deqn{p | X \sim Beta(\alpha + X, \beta + N - X)}
#' Posterior mean:
#' \deqn{E[p|X] = (\alpha + X) / (\alpha + \beta + N)}
#'
#' This provides principled shrinkage of extreme observed proportions (0/1)
#' toward \eqn{\pi}, with strength implied by \eqn{\rho}.
#'
#' @param count numeric vector of counts X (0..N)
#' @param library_size numeric vector of N (>0)
#' @param pi numeric vector of prior mean in [0,1]
#' @param rho numeric vector of ICC in (0,1)
#' @param eps clamp for pi and rho away from boundaries
#'
#' @return numeric vector posterior mean in (0,1)
bb_posterior_mean_p <- function(count, library_size, pi, rho, eps = 1e-8) {
  if (!is.numeric(count)) stop("count must be numeric.")
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(pi)) stop("pi must be numeric.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")

  # Clamp away from exact boundaries for numerical stability
  rho <- pmax(pmin(rho, 1 - eps), eps)
  pi <- .clamp01(pi, eps = eps)

  # Guard counts within [0, N]
  count <- pmax(0, pmin(count, library_size))

  phi <- (1 / rho) - 1
  alpha <- pi * phi
  beta <- (1 - pi) * phi

  (alpha + count) / (alpha + beta + library_size)
}

#' Posterior mean of arcsin-sqrt transformed proportion under Beta posterior (delta approx)
#'
#' Instead of shrinking on the proportion scale (p), we approximate shrinkage directly on
#' the arcsin-sqrt scale by approximating \(E[g(p)\mid X]\) where \(g(p)=\arcsin\sqrt{p}\).
#' This introduces explicit **count-dependence** (via the Beta posterior) so that
#' low-information observations (e.g. very small X, or boundary outcomes X∈{0,N})
#' are regularized more strongly.
#'
#' Approximation (second-order delta method around posterior mean m):
#' \deqn{E[g(p)\mid X] \approx g(m) + \tfrac{1}{2} g''(m)\,Var(p\mid X)}
#'
#' @param count numeric vector of counts X (0..N)
#' @param library_size numeric vector of N (>0)
#' @param pi numeric vector prior mean in [0,1]
#' @param rho numeric vector ICC in (0,1)
#' @param eps clamp for numerical stability
#'
#' @return numeric vector approx \(E[\arcsin(\sqrt{p})\mid X]\)
bb_posterior_mean_arcsin_sqrt <- function(count, library_size, pi, rho, eps = 1e-8) {
  if (!is.numeric(count)) stop("count must be numeric.")
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(pi)) stop("pi must be numeric.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")

  rho <- pmax(pmin(rho, 1 - eps), eps)
  pi <- .clamp01(pi, eps = eps)
  count <- pmax(0, pmin(count, library_size))

  phi <- (1 / rho) - 1
  alpha <- pi * phi
  beta <- (1 - pi) * phi

  a_post <- alpha + count
  b_post <- beta + (library_size - count)
  ab_post <- a_post + b_post

  m <- a_post / ab_post
  m <- .clamp01(m, eps = eps)

  # Var(p | X) for Beta(a,b): ab / ((a+b)^2 (a+b+1))
  v <- (a_post * b_post) / ((ab_post^2) * (ab_post + 1))

  g <- asin(sqrt(m))
  # g''(p) = (2p-1) / (4*(p(1-p))^(3/2))
  gpp <- (2 * m - 1) / (4 * (m * (1 - m))^(3/2))

  g + 0.5 * gpp * v
}

#' Approximate SD on arcsin-sqrt scale (pi-aware; second-order delta method)
#'
#' First-order delta method cancels pi dependence:
#' \deqn{Var(asin(sqrt(p_hat))) \approx (1 + (N-1)rho) / (4N)}
#'
#' Near boundaries (pi close to 0/1), first-order can be inaccurate. This
#' pi-aware approximation adds a second-order term:
#' \deqn{Var(g(p_hat)) \approx g'(pi)^2 Var(p_hat) + 0.5 g''(pi)^2 Var(p_hat)^2}
#'
#' where g(p) = asin(sqrt(p)).
#'
#' @param pi numeric vector in [0,1]
#' @param library_size numeric vector of N (>0)
#' @param rho numeric vector in (0,1)
#' @param eps clamp for pi away from 0/1
#'
#' @return numeric vector SD approximation
bb_arcsin_sqrt_sd_piaware <- function(pi, library_size, rho, eps = 1e-8) {
  pi <- .clamp01(pi, eps = eps)
  v_p <- bb_var_p_hat(pi = pi, library_size = library_size, rho = rho, eps = eps)

  # g'(p) = 1/(2*sqrt(p(1-p)))
  gp <- 1 / (2 * sqrt(pi * (1 - pi)))
  # g''(p) = (2p-1) / (4*(p(1-p))^(3/2))
  gpp <- (2 * pi - 1) / (4 * (pi * (1 - pi))^(3/2))

  v1 <- (gp^2) * v_p
  v2 <- 0.5 * (gpp^2) * (v_p^2)
  sqrt(v1 + v2)
}

#' Association-adjust residuals in arcsin-sqrt space using beta-binomial SD
#'
#' This is the "correct" variance-scale adjustment when residuals are computed as:
#' \deqn{r_{ij} = \arcsin\sqrt{\hat p_{ij}} - \arcsin\sqrt{\pi_{gj}}}
#'
#' We assume a mean–dispersion association in unconstrained space:
#' \deqn{\log(\sigma_{ij}) = -k \cdot \mu_{gj} + c}
#'
#' and adjust residuals by the *relative* arcsin-sqrt SD implied by the association:
#' \deqn{r^*_{ij} = r_{ij} \cdot \frac{\mathrm{SD}_0(N_i)}{\mathrm{SD}(N_i, \rho_{ij})}}
#'
#' where \eqn{\mathrm{SD}_0(N_i)} is the SD at \eqn{\mu = 0} (i.e. \eqn{\sigma_0 = \exp(c)}),
#' so this is "rotation around the origin" (no global inflation/deflation).
#'
#' @param residual numeric residuals on arcsin-sqrt scale
#' @param mu_inv_softmax numeric unconstrained predictors (same shape as residual)
#' @param log_dispersion_assoc numeric scalar k
#' @param library_size numeric vector of N (same shape as residual, or recycled)
#' @param log_sigma_intercept optional numeric scalar c. If NULL, it will be estimated from
#'   `log_sigma_observed` as mean(log_sigma_observed + k * mu_inv_softmax).
#' @param log_sigma_observed optional numeric vector of observed log_sigma (same shape as residual).
#'   Used only if log_sigma_intercept is NULL.
#'
#' @return numeric adjusted residuals
adjust_residual_assoc_arcsin_bb <- function(
  residual,
  mu_inv_softmax,
  log_dispersion_assoc,
  library_size,
  log_sigma_intercept
) {
  if (!is.numeric(residual)) stop("residual must be numeric.")
  if (!is.numeric(mu_inv_softmax)) stop("mu_inv_softmax must be numeric.")
  if (!is.numeric(log_dispersion_assoc) || length(log_dispersion_assoc) != 1) {
    stop("log_dispersion_assoc must be a numeric scalar.")
  }
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")

  # Model: log_sigma = -k * mu + c
  # Here we assume `log_sigma_intercept` (c) is KNOWN (e.g. from simulation inputs).
  # This keeps the adjustment fully algebraic (no estimation).

  # Baseline at mu = 0
  sigma0 <- exp(log_sigma_intercept)
  rho0 <- sigma_to_rho(sigma0)

  # Expected sigma at each mu under the model
  log_sigma_exp <- -log_dispersion_assoc * mu_inv_softmax + log_sigma_intercept
  sigma_exp <- exp(log_sigma_exp)
  rho_exp <- sigma_to_rho(sigma_exp)

  var0 <- (1 + (library_size - 1) * rho0) / (4 * library_size)
  var_exp <- (1 + (library_size - 1) * rho_exp) / (4 * library_size)

  residual * sqrt(var0 / var_exp)
}

#' "inv_softmax" transform: map unconstrained log-linear predictors to proportions
#'
#' In this repo's naming, "mu_inv_softmax" refers to an *unconstrained* vector
#' (log-linear predictors). Applying "inv_softmax" converts it to proportions.
#' Mathematically, this is the SOFTMAX transform.
#'
#' @param x numeric vector or matrix of unconstrained values (real-valued).
#'   - If a vector: returns a vector of the same length that sums to 1.
#'   - If a matrix: applies row-wise; each row sums to 1.
#' @param eps numeric scalar. If > 0, probabilities are clamped to [eps, 1-eps].
#'
#' @return numeric vector or matrix of probabilities in (0,1), same shape as x.
inv_softmax <- function(x, eps = 0) {
  if (!is.numeric(x)) stop("x must be numeric.")

  softmax_vec <- function(v) {
    if (anyNA(v)) stop("x contains NA; please handle missing values before softmax.")
    v <- v - max(v) # numerical stability
    ex <- exp(v)
    out <- ex / sum(ex)
    if (eps > 0) out <- .clamp01(out, eps = eps)
    out
  }

  if (is.matrix(x)) {
    out <- t(apply(x, 1, softmax_vec))
    dimnames(out) <- dimnames(x)
    return(out)
  }

  if (is.vector(x)) {
    return(softmax_vec(as.numeric(x)))
  }

  stop("x must be a numeric vector or matrix.")
}


