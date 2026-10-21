// Fisher scoring (IRLS) for GLMs, RcppArmadillo engine.
//
// This mirrors fs_fit_r() in R/fit_r.R step for step, so the two engines
// give the same answer to machine precision. Family and link are passed as
// integer codes:
//   family: 1 = binomial, 2 = poisson, 3 = gaussian
//   link:   1 = logit, 2 = probit, 3 = cloglog, 4 = log, 5 = identity, 6 = sqrt

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <cfloat>

using namespace Rcpp;

namespace {

const double EPS = DBL_EPSILON;

// ---- inverse link g^{-1}(eta) ------------------------------------------
double linkinv(double eta, int link) {
  switch (link) {
  case 1: {  // logit, thresholds as in stats:::C_logit_linkinv
    double tmp = (eta < -30.0) ? EPS : ((eta > 30.0) ? 1.0 / EPS : std::exp(eta));
    return tmp / (1.0 + tmp);
  }
  case 2: {  // probit
    double thresh = -R::qnorm(EPS, 0.0, 1.0, 1, 0);
    double e = std::min(std::max(eta, -thresh), thresh);
    return R::pnorm(e, 0.0, 1.0, 1, 0);
  }
  case 3:  // cloglog
    return std::max(std::min(-std::expm1(-std::exp(eta)), 1.0 - EPS), EPS);
  case 4:  // log
    return std::max(std::exp(eta), EPS);
  case 5:  // identity
    return eta;
  case 6:  // sqrt
    return eta * eta;
  }
  stop("unknown link code");
}

// ---- d mu / d eta --------------------------------------------------------
double mu_eta(double eta, int link) {
  switch (link) {
  case 1: {
    if (std::fabs(eta) > 30.0) return EPS;
    double ex = std::exp(eta);
    return ex / ((1.0 + ex) * (1.0 + ex));
  }
  case 2:
    return std::max(R::dnorm(eta, 0.0, 1.0, 0), EPS);
  case 3: {
    double e = std::min(eta, 700.0);
    return std::max(std::exp(e) * std::exp(-std::exp(e)), EPS);
  }
  case 4:
    return std::max(std::exp(eta), EPS);
  case 5:
    return 1.0;
  case 6:
    return 2.0 * eta;
  }
  stop("unknown link code");
}

// ---- variance function V(mu) ---------------------------------------------
double variance(double mu, int family) {
  switch (family) {
  case 1: return mu * (1.0 - mu);
  case 2: return mu;
  case 3: return 1.0;
  }
  stop("unknown family code");
}

inline double ylogy(double y, double mu) { return (y > 0.0) ? y * std::log(y / mu) : 0.0; }

// ---- total deviance -------------------------------------------------------
double deviance(const arma::vec& y, const arma::vec& mu, const arma::vec& wt, int family) {
  double dev = 0.0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    switch (family) {
    case 1: dev += 2.0 * wt[i] * (ylogy(y[i], mu[i]) + ylogy(1.0 - y[i], 1.0 - mu[i])); break;
    case 2: dev += 2.0 * wt[i] * (ylogy(y[i], mu[i]) - (y[i] - mu[i])); break;
    case 3: dev += wt[i] * (y[i] - mu[i]) * (y[i] - mu[i]); break;
    }
  }
  return dev;
}

bool valid(const arma::vec& eta, const arma::vec& mu, int family, int link) {
  for (arma::uword i = 0; i < eta.n_elem; ++i) {
    if (link == 6 && !(std::isfinite(eta[i]) && eta[i] > 0.0)) return false;
    if (family == 1 && !(std::isfinite(mu[i]) && mu[i] > 0.0 && mu[i] < 1.0)) return false;
    if (family == 2 && !(std::isfinite(mu[i]) && mu[i] > 0.0)) return false;
  }
  return true;
}

arma::vec linkinv_vec(const arma::vec& eta, int link) {
  arma::vec mu(eta.n_elem);
  for (arma::uword i = 0; i < eta.n_elem; ++i) mu[i] = linkinv(eta[i], link);
  return mu;
}

}  // namespace

//' Fisher scoring in C++
//'
//' Low-level RcppArmadillo engine with the same algorithm and arguments as
//' [fs_fit_r()], except that the family and link are passed as integer codes.
//' Most users should call [fs_glm()] with `engine = "cpp"`.
//'
//' @param X Model matrix (full column rank).
//' @param y Response (proportions for binomial).
//' @param weights Prior weights.
//' @param offset Offset vector.
//' @param eta_start Starting linear predictor.
//' @param start Starting coefficients (length 0 if none).
//' @param family_code,link_code Integer codes, see `fs_family()`.
//' @param epsilon,maxit,trace Control settings, see [fs_control()].
//' @return A list with the same components as [fs_fit_r()].
//' @keywords internal
//' @export
// [[Rcpp::export]]
List fs_fit_cpp(const arma::mat& X, const arma::vec& y, const arma::vec& weights,
                const arma::vec& offset, const arma::vec& eta_start,
                const arma::vec& start, int family_code, int link_code,
                double epsilon, int maxit, bool trace) {
  const arma::uword n = X.n_rows, p = X.n_cols;

  arma::vec eta = eta_start;
  arma::vec mu = linkinv_vec(eta, link_code);
  if (!valid(eta, mu, family_code, link_code))
    stop("cannot find valid starting values: please specify some");
  double devold = deviance(y, mu, weights, family_code);
  double dev = devold;

  bool have_coefold = start.n_elem == p;
  arma::vec coefold = have_coefold ? start : arma::vec(p, arma::fill::zeros);
  arma::vec beta(p), coef(p);
  coef.fill(NA_REAL);
  arma::vec mu_eta_v(n), z(n), sw(n);
  arma::mat Xw(n, p), R;
  arma::rowvec colnorm(p);
  bool conv = false, boundary = false;
  int iter = 0;
  std::vector<double> trace_dev;
  std::vector<arma::vec> trace_coef;

  for (iter = 1; iter <= maxit; ++iter) {
    R_CheckUserInterrupt();

    // 1. Working weights and working response
    for (arma::uword i = 0; i < n; ++i) {
      mu_eta_v[i] = mu_eta(eta[i], link_code);
      z[i] = (eta[i] - offset[i]) + (y[i] - mu[i]) / mu_eta_v[i];
      sw[i] = std::sqrt(weights[i] * mu_eta_v[i] * mu_eta_v[i] / variance(mu[i], family_code));
    }

    // 2. Weighted least squares by QR. Columns of W^{1/2} X are first scaled to
    // unit norm (W^{1/2} X = Xs D), so the fit does not depend on predictor
    // units: Xs = Q R, R gamma = Q' W^{1/2} z, beta = D^{-1} gamma.
    Xw = X.each_col() % sw;
    colnorm = arma::sqrt(arma::sum(arma::square(Xw), 0));
    if (arma::any(colnorm == 0.0)) stop("X'WX is singular at iteration %d", iter);
    Xw.each_row() /= colnorm;
    arma::mat Q;
    if (!arma::qr_econ(Q, R, Xw)) stop("QR decomposition failed at iteration %d", iter);
    // Column j is (numerically) dependent on the earlier ones when |R_jj|, the
    // norm of the part of its unit-norm column orthogonal to them, is tiny.
    // Same tolerance as fs_fit_r().
    for (arma::uword j = 0; j < p; ++j)
      if (!(std::fabs(R(j, j)) > std::min(1e-7, epsilon / 1000.0)))
        stop("X'WX is singular at iteration %d", iter);
    beta = arma::solve(arma::trimatu(R), Q.t() * (z % sw)) / colnorm.t();

    // 3. New linear predictor, mean and deviance
    eta = X * beta + offset;
    mu = linkinv_vec(eta, link_code);
    dev = deviance(y, mu, weights, family_code);
    if (trace) Rprintf("Deviance = %.10g Iterations - %d\n", dev, iter);

    // 4. Step-halving for a non-finite deviance or invalid eta/mu
    if (!std::isfinite(dev)) {
      if (!have_coefold) stop("no valid set of coefficients has been found: please supply starting values");
      int ii = 1;
      while (!std::isfinite(dev)) {
        if (ii > maxit) stop("inner loop 1; cannot correct step size");
        ++ii;
        beta = (beta + coefold) / 2.0;
        eta = X * beta + offset;
        mu = linkinv_vec(eta, link_code);
        dev = deviance(y, mu, weights, family_code);
      }
      boundary = true;
    }
    if (!valid(eta, mu, family_code, link_code)) {
      if (!have_coefold) stop("no valid set of coefficients has been found: please supply starting values");
      int ii = 1;
      while (!valid(eta, mu, family_code, link_code)) {
        if (ii > maxit) stop("inner loop 2; cannot correct step size");
        ++ii;
        beta = (beta + coefold) / 2.0;
        eta = X * beta + offset;
        mu = linkinv_vec(eta, link_code);
      }
      boundary = true;
      dev = deviance(y, mu, weights, family_code);
    }

    trace_dev.push_back(dev);
    trace_coef.push_back(beta);

    // 5. Convergence check (same rule as glm)
    if (std::fabs(dev - devold) / (std::fabs(dev) + 0.1) < epsilon) {
      conv = true;
      coef = beta;
      break;
    }
    devold = dev;
    coef = coefold = beta;
    have_coefold = true;
  }
  if (iter > maxit) iter = maxit;

  // Unscaled covariance (X'WX)^{-1} = D^{-1} R^{-1} R^{-T} D^{-1} from the
  // final WLS step
  arma::mat Rinv = arma::inv(arma::trimatu(R));
  Rinv.each_col() /= colnorm.t();
  arma::mat cov_unscaled = Rinv * Rinv.t();

  arma::mat tc(trace_coef.size(), p);
  for (size_t k = 0; k < trace_coef.size(); ++k) tc.row(k) = trace_coef[k].t();

  return List::create(
    _["coefficients"] = NumericVector(coef.begin(), coef.end()),
    _["linear.predictors"] = NumericVector(eta.begin(), eta.end()),
    _["fitted.values"] = NumericVector(mu.begin(), mu.end()),
    _["deviance"] = dev,
    _["iter"] = iter,
    _["converged"] = conv,
    _["boundary"] = boundary,
    _["working.weights"] = NumericVector(sw.begin(), sw.end()) * NumericVector(sw.begin(), sw.end()),
    _["cov.unscaled"] = cov_unscaled,
    _["trace"] = List::create(_["deviance"] = wrap(trace_dev), _["coefficients"] = tc));
}
