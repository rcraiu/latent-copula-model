// Fused C++ kernel for the gamma-Gaussian latent variable copula model.
//
// Profiling the pure-R sampler showed that ~82% of the time went into three
// places: the copula log density (46%, of which clamp01's pmin/pmax alone was
// 21%), the gamma margin (29%), and the normal margin (7%). All three are
// per-observation and were being evaluated through separate vectorised R calls,
// each allocating several length-n temporaries. Fusing them into one loop
// removes the allocations and the repeated passes over memory.
//
// One algebraic saving is also taken. For the Gaussian copula the kernel needs
// z2 = Phi^{-1}(F_y(y|w)), and F_y is itself normal, so the pnorm/qnorm round
// trip is the identity: z2 = (y - mu_y - lam_y w)/sigma directly. The bound
// ZCLAMP reproduces the clamping that the R code applies on the u scale.

#include <Rcpp.h>
using namespace Rcpp;

static const double CLAMP  = 1e-12;
static const double ZCLAMP = 7.0339021663017056;  // qnorm(1 - 1e-12)

enum Family { GAUSSIAN = 0, CLAYTON = 1, FRANK = 2, GUMBEL = 3, INDEP = 4 };
enum What   { JOINT = 0, Q_GIVEN_Y = 1, Y_GIVEN_Q = 2 };

static inline double clamp01(double u) {
  if (u < CLAMP) return CLAMP;
  if (u > 1.0 - CLAMP) return 1.0 - CLAMP;
  return u;
}

static inline double theta_of(double w, double xi0, double xi1,
                              int family, double gscale) {
  double lin = gscale * (xi0 + xi1 * w);
  switch (family) {
    case GAUSSIAN: {
      double t = std::tanh(lin);
      if (t >  1.0 - 1e-8) t =  1.0 - 1e-8;
      if (t < -1.0 + 1e-8) t = -1.0 + 1e-8;
      return t;
    }
    case FRANK:   return (std::fabs(lin) < 1e-6) ? 1e-6 : lin;
    case GUMBEL:  return 1.0 + std::exp(lin);
    case CLAYTON: return std::exp(lin);
    default:      return 0.0;
  }
}

// log copula density. For the Gaussian family z2 is supplied directly and u2 is
// not needed; the other families work on the u scale.
static inline double log_cop(double u1, double u2, double z1, double z2,
                             double th, int family) {
  switch (family) {
    case GAUSSIAN: {
      double r = th;
      if (r >  1.0 - 1e-8) r =  1.0 - 1e-8;
      if (r < -1.0 + 1e-8) r = -1.0 + 1e-8;
      double om = 1.0 - r * r;
      return -0.5 * std::log(om)
             - (z1 * z1 - 2.0 * r * z1 * z2 + z2 * z2) / (2.0 * om)
             + 0.5 * (z1 * z1 + z2 * z2);
    }
    case CLAYTON: {
      double t = th < 1e-10 ? 1e-10 : th;
      double lu1 = std::log(u1), lu2 = std::log(u2);
      double A = std::pow(u1, -t) + std::pow(u2, -t) - 1.0;
      if (!(A > 0.0)) return R_NegInf;
      return std::log1p(t) - (1.0 + t) * (lu1 + lu2) - (2.0 + 1.0 / t) * std::log(A);
    }
    case GUMBEL: {
      double t = th < 1.0 + 1e-10 ? 1.0 + 1e-10 : th;
      double x = -std::log(u1), yy = -std::log(u2);
      double S = std::pow(x, t) + std::pow(yy, t);
      double A = std::pow(S, 1.0 / t);
      return -A - std::log(u1) - std::log(u2)
             + (t - 1.0) * (std::log(x) + std::log(yy))
             + (1.0 - 2.0 * t) * std::log(A) + std::log(A + t - 1.0);
    }
    case FRANK: {
      double t = std::fabs(th) < 1e-8 ? 1e-8 : th;
      double g1  = std::expm1(-t);
      double num = -t * g1 * (1.0 + std::expm1(-t * (u1 + u2)));
      double den = std::expm1(-t * u1) * std::expm1(-t * u2) + g1;
      if (!(num > 0.0)) return R_NegInf;
      return std::log(num) - 2.0 * std::log(std::fabs(den));
    }
    default: return 0.0;
  }
}

//' Per-observation log density of the surrogate pair given the latent variable.
//'
//' @param what 0 = joint f(q,y|w), 1 = f(q|y,w), 2 = f(y|q,w)
// [[Rcpp::export]]
NumericVector lvc_logdens_cpp(const NumericVector& q, const NumericVector& y,
                              const NumericVector& W,
                              double muq, double lamq, double alpha,
                              double muy, double lamy, double sigma,
                              double xi0, double xi1,
                              int family, double gscale, int what) {
  const int n = q.size();
  NumericVector out(n);
  const double lsig = std::log(sigma);
  const bool need_u2 = (family != GAUSSIAN) && (family != INDEP);

  for (int i = 0; i < n; ++i) {
    const double w = W[i];

    // gamma margin, mean exp(muq + lamq w), shape alpha
    const double scale = std::exp(muq + lamq * w) / alpha;
    const double lf1 = R::dgamma(q[i], alpha, scale, 1);
    double u1 = clamp01(R::pgamma(q[i], alpha, scale, 1, 0));

    // normal margin
    double z2 = (y[i] - muy - lamy * w) / sigma;
    const double lf2 = R::dnorm(z2, 0.0, 1.0, 1) - lsig;
    if (z2 >  ZCLAMP) z2 =  ZCLAMP;
    if (z2 < -ZCLAMP) z2 = -ZCLAMP;

    double lc = 0.0;
    if (family != INDEP) {
      const double th = theta_of(w, xi0, xi1, family, gscale);
      const double z1 = R::qnorm(u1, 0.0, 1.0, 1, 0);
      const double u2 = need_u2 ? clamp01(R::pnorm(z2, 0.0, 1.0, 1, 0)) : 0.0;
      lc = log_cop(u1, u2, z1, z2, th, family);
    }

    double v;
    switch (what) {
      case Q_GIVEN_Y: v = lf1 + lc; break;
      case Y_GIVEN_Q: v = lf2 + lc; break;
      default:        v = lf1 + lf2 + lc; break;
    }
    out[i] = R_finite(v) ? v : R_NegInf;
  }
  return out;
}

//' Sum of the per-observation log joint density, without allocating the vector.
// [[Rcpp::export]]
double lvc_loglik_cpp(const NumericVector& q, const NumericVector& y,
                      const NumericVector& W,
                      double muq, double lamq, double alpha,
                      double muy, double lamy, double sigma,
                      double xi0, double xi1,
                      int family, double gscale) {
  const int n = q.size();
  const double lsig = std::log(sigma);
  const bool need_u2 = (family != GAUSSIAN) && (family != INDEP);
  double acc = 0.0;

  for (int i = 0; i < n; ++i) {
    const double w = W[i];
    const double scale = std::exp(muq + lamq * w) / alpha;
    const double lf1 = R::dgamma(q[i], alpha, scale, 1);
    double u1 = clamp01(R::pgamma(q[i], alpha, scale, 1, 0));

    double z2 = (y[i] - muy - lamy * w) / sigma;
    const double lf2 = R::dnorm(z2, 0.0, 1.0, 1) - lsig;
    if (z2 >  ZCLAMP) z2 =  ZCLAMP;
    if (z2 < -ZCLAMP) z2 = -ZCLAMP;

    double lc = 0.0;
    if (family != INDEP) {
      const double th = theta_of(w, xi0, xi1, family, gscale);
      const double z1 = R::qnorm(u1, 0.0, 1.0, 1, 0);
      const double u2 = need_u2 ? clamp01(R::pnorm(z2, 0.0, 1.0, 1, 0)) : 0.0;
      lc = log_cop(u1, u2, z1, z2, th, family);
    }

    const double v = lf1 + lf2 + lc;
    if (!R_finite(v)) return R_NegInf;
    acc += v;
  }
  return acc;
}

//' All three predictive densities in one pass.
//'
//' The criteria need f(q,y|w), f(q|y,w) and f(y|q,w) at every retained
//' iteration. Evaluating them through three separate calls recomputes the
//' copula density three times; here the shared pieces are computed once.
//' Returns an n x 3 matrix with columns (joint, q|y, y|q).
// [[Rcpp::export]]
NumericMatrix lvc_logdens3_cpp(const NumericVector& q, const NumericVector& y,
                               const NumericVector& W,
                               double muq, double lamq, double alpha,
                               double muy, double lamy, double sigma,
                               double xi0, double xi1,
                               int family, double gscale) {
  const int n = q.size();
  NumericMatrix out(n, 3);
  const double lsig = std::log(sigma);
  const bool need_u2 = (family != GAUSSIAN) && (family != INDEP);

  for (int i = 0; i < n; ++i) {
    const double w = W[i];
    const double scale = std::exp(muq + lamq * w) / alpha;
    const double lf1 = R::dgamma(q[i], alpha, scale, 1);
    double u1 = clamp01(R::pgamma(q[i], alpha, scale, 1, 0));

    double z2 = (y[i] - muy - lamy * w) / sigma;
    const double lf2 = R::dnorm(z2, 0.0, 1.0, 1) - lsig;
    if (z2 >  ZCLAMP) z2 =  ZCLAMP;
    if (z2 < -ZCLAMP) z2 = -ZCLAMP;

    double lc = 0.0;
    if (family != INDEP) {
      const double th = theta_of(w, xi0, xi1, family, gscale);
      const double z1 = R::qnorm(u1, 0.0, 1.0, 1, 0);
      const double u2 = need_u2 ? clamp01(R::pnorm(z2, 0.0, 1.0, 1, 0)) : 0.0;
      lc = log_cop(u1, u2, z1, z2, th, family);
    }

    const double j = lf1 + lf2 + lc, a = lf1 + lc, b = lf2 + lc;
    out(i, 0) = R_finite(j) ? j : R_NegInf;
    out(i, 1) = R_finite(a) ? a : R_NegInf;
    out(i, 2) = R_finite(b) ? b : R_NegInf;
  }
  return out;
}
