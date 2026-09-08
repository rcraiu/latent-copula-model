// Full Gibbs sampler for the gamma-Gaussian latent variable copula model.
//
// Why the whole loop and not just a density kernel: timing the primitives at
// n = 3948 (300 passes) gives pgamma 0.250s, dgamma 0.098s, qnorm 0.025s,
// pnorm 0.031s, clamp 0.009s. The gamma margin is roughly 70% of the cost and
// is already R's own compiled code, so a fused density kernel only buys ~1.3x.
// The real savings are structural, and they need the cache to live inside the
// loop:
//
//   * The gamma pieces (lf1, u1, z1, lu1) depend only on (muq, lamq, alpha, W),
//     so they are recomputed only when that block or W moves -- never for the
//     y-margin or copula proposals.
//   * z1 = Phi^{-1}(u1) and the logs of the pseudo-observations are cached
//     rather than recomputed inside every copula evaluation. The pure-R version
//     called qnorm twice per copula density, 14 times per iteration.
//   * For the Gaussian copula the normal margin gives z2 = (y - mu_y - lam_y w)/sigma
//     analytically, so its pnorm/qnorm round trip disappears entirely.
//   * The three predictive densities needed by the criteria share one pass.
//
// The sampler otherwise reproduces the R implementation step for step:
// Andrieu & Thoms (2008) Algorithm 6 (componentwise adaptive scaling) for the
// parameter blocks and Algorithm 4 (global adaptive scaling) for each W_i, both
// targeting an acceptance rate of 0.44.

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

static const double CLAMP  = 1e-12;
static const double ZCLAMP = 7.0339021663017056;   // qnorm(1 - 1e-12)
static const double LOG_SQRT_2PI = 0.9189385332046727;

enum Family { GAUSSIAN = 0, CLAYTON = 1, FRANK = 2, GUMBEL = 3, INDEP = 4 };

// Optimal acceptance rates for a random walk in d dimensions (Roberts, Gelman &
// Gilks 1997): 0.44 for a scalar update, falling to 0.234 as d grows. Indexed by
// block length; only lengths 2 and 3 occur here.
static const double BLOCK_TARGET[4] = {0.44, 0.44, 0.352, 0.316};

static inline double clamp01(double u) {
  if (u < CLAMP) return CLAMP;
  if (u > 1.0 - CLAMP) return 1.0 - CLAMP;
  return u;
}

static inline double theta_of(double w, double xi0, double xi1,
                              int family, double gscale) {
  const double lin = gscale * (xi0 + xi1 * w);
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

// Log copula density from cached pieces: z1, z2 are the normal scores and
// lu1, lu2 the logs of the pseudo-observations.
static inline double log_cop(double u1, double u2, double lu1, double lu2,
                             double z1, double z2, double th, int family) {
  switch (family) {
    case GAUSSIAN: {
      double r = th;
      if (r >  1.0 - 1e-8) r =  1.0 - 1e-8;
      if (r < -1.0 + 1e-8) r = -1.0 + 1e-8;
      const double om = 1.0 - r * r;
      return -0.5 * std::log(om)
             - (z1 * z1 - 2.0 * r * z1 * z2 + z2 * z2) / (2.0 * om)
             + 0.5 * (z1 * z1 + z2 * z2);
    }
    case CLAYTON: {
      const double t = th < 1e-10 ? 1e-10 : th;
      const double A = std::exp(-t * lu1) + std::exp(-t * lu2) - 1.0;
      if (!(A > 0.0)) return R_NegInf;
      return std::log1p(t) - (1.0 + t) * (lu1 + lu2) - (2.0 + 1.0 / t) * std::log(A);
    }
    case GUMBEL: {
      const double t = th < 1.0 + 1e-10 ? 1.0 + 1e-10 : th;
      const double x = -lu1, yy = -lu2;
      const double S = std::pow(x, t) + std::pow(yy, t);
      const double A = std::pow(S, 1.0 / t);
      return -A - lu1 - lu2 + (t - 1.0) * (std::log(x) + std::log(yy))
             + (1.0 - 2.0 * t) * std::log(A) + std::log(A + t - 1.0);
    }
    case FRANK: {
      const double t = std::fabs(th) < 1e-8 ? 1e-8 : th;
      const double g1  = std::expm1(-t);
      const double num = -t * g1 * (1.0 + std::expm1(-t * (u1 + u2)));
      const double den = std::expm1(-t * u1) * std::expm1(-t * u2) + g1;
      if (!(num > 0.0)) return R_NegInf;
      return std::log(num) - 2.0 * std::log(std::fabs(den));
    }
    default: return 0.0;
  }
}

// Cached per-observation pieces of the joint density.
struct Cache {
  std::vector<double> lf1, u1, lu1, z1;   // gamma margin
  std::vector<double> lf2, u2, lu2, z2;   // normal margin
  std::vector<double> th;                 // copula parameter
  Cache(int n) : lf1(n), u1(n), lu1(n), z1(n),
                 lf2(n), u2(n), lu2(n), z2(n), th(n) {}
};

struct Model {
  const double *q, *y;
  int n, family;
  double gscale;
  bool need_u;                            // do the copula formulas need u/log u?

  void fill_q(std::vector<double>& lf1, std::vector<double>& u1,
              std::vector<double>& lu1, std::vector<double>& z1,
              const double* W, double muq, double lamq, double alpha) const {
    for (int i = 0; i < n; ++i) {
      const double scale = std::exp(muq + lamq * W[i]) / alpha;
      lf1[i] = R::dgamma(q[i], alpha, scale, 1);
      const double u = clamp01(R::pgamma(q[i], alpha, scale, 1, 0));
      u1[i]  = u;
      if (family == GAUSSIAN) z1[i] = R::qnorm(u, 0.0, 1.0, 1, 0);
      if (need_u)             lu1[i] = std::log(u);
    }
  }

  void fill_y(std::vector<double>& lf2, std::vector<double>& u2,
              std::vector<double>& lu2, std::vector<double>& z2,
              const double* W, double muy, double lamy, double sigma) const {
    const double lsig = std::log(sigma);
    for (int i = 0; i < n; ++i) {
      double z = (y[i] - muy - lamy * W[i]) / sigma;
      lf2[i] = -0.5 * z * z - LOG_SQRT_2PI - lsig;
      if (z >  ZCLAMP) z =  ZCLAMP;
      if (z < -ZCLAMP) z = -ZCLAMP;
      z2[i] = z;                          // Phi^{-1}(F_y) is the identity here
      if (need_u) {
        const double u = clamp01(R::pnorm(z, 0.0, 1.0, 1, 0));
        u2[i]  = u;
        lu2[i] = std::log(u);
      }
    }
  }

  void fill_th(std::vector<double>& th, const double* W,
               double xi0, double xi1) const {
    if (family == INDEP) { std::fill(th.begin(), th.end(), 0.0); return; }
    for (int i = 0; i < n; ++i) th[i] = theta_of(W[i], xi0, xi1, family, gscale);
  }

  inline double ld(const Cache& C, int i) const {
    const double lc = (family == INDEP) ? 0.0
      : log_cop(C.u1[i], C.u2[i], C.lu1[i], C.lu2[i], C.z1[i], C.z2[i],
                C.th[i], family);
    const double v = C.lf1[i] + C.lf2[i] + lc;
    return R_finite(v) ? v : R_NegInf;
  }

  double loglik(const Cache& C) const {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
      const double v = ld(C, i);
      if (!R_finite(v)) return R_NegInf;
      s += v;
    }
    return s;
  }
};

// ---------------------------------------------------------------- prior ----
static double logprior(const double* p, bool has_xi, double s2_mu, double s2_lam,
                       double s2_xi, double a_alpha, double b_alpha,
                       double a_sig, double b_sig) {
  // p = (muq, lamq, lalpha, muy, lamy, lsigma, xi0, xi1)
  double lp = R::dnorm(p[0], 0.0, std::sqrt(s2_mu),  1)
            + R::dnorm(p[1], 0.0, std::sqrt(s2_lam), 1)
            + R::dnorm(p[3], 0.0, std::sqrt(s2_mu),  1)
            + R::dnorm(p[4], 0.0, std::sqrt(s2_lam), 1)
            + R::dgamma(std::exp(p[2]), a_alpha, 1.0 / b_alpha, 1) + p[2]
            + (-a_sig - 1.0) * (2.0 * p[5]) - b_sig / std::exp(2.0 * p[5])
            + 2.0 * p[5];
  if (has_xi) lp += R::dnorm(p[6], 0.0, std::sqrt(s2_xi), 1)
                  + R::dnorm(p[7], 0.0, std::sqrt(s2_xi), 1);
  return lp;
}

//' Gibbs sampler for the gamma-Gaussian latent variable copula model.
//'
//' @param blockwise if true the three marginal/copula parameters of a block are
//'   proposed jointly from the adapted covariance (one gamma evaluation per
//'   block per iteration) instead of one component at a time.
// [[Rcpp::export]]
List run_lvc_cpp(NumericVector q, NumericVector y, NumericMatrix X,
                 int family, int M, int burn, double gscale,
                 NumericVector p_init, NumericVector beta_init,
                 NumericMatrix betapost_cov, NumericMatrix betapost_cholT,
                 List prior, int thin_W, bool track_criteria, int verbose,
                 bool blockwise = false, bool cross_block = false) {
  const int n = q.size(), pd = X.ncol();
  const bool has_xi = (family != INDEP);
  const double TARGET = 0.44;

  Model mod;
  mod.q = &q[0]; mod.y = &y[0]; mod.n = n; mod.family = family;
  mod.gscale = gscale;
  mod.need_u = (family == CLAYTON || family == FRANK || family == GUMBEL);

  const double s2_mu = prior["s2_mu"], s2_lam = prior["s2_lam"],
               s2_xi = prior["s2_xi"], a_alpha = prior["a_alpha"],
               b_alpha = prior["b_alpha"], a_sig = prior["a_sig"],
               b_sig = prior["b_sig"];

  std::vector<double> p(p_init.begin(), p_init.end());
  std::vector<double> beta(beta_init.begin(), beta_init.end());
  std::vector<double> W(n);
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = 0; k < pd; ++k) s += X(i, k) * beta[k];
    W[i] = s;
  }

  Cache C(n), S(n);                       // current cache and scratch
  mod.fill_q(C.lf1, C.u1, C.lu1, C.z1, W.data(), p[0], p[1], std::exp(p[2]));
  mod.fill_y(C.lf2, C.u2, C.lu2, C.z2, W.data(), p[3], p[4], std::exp(p[5]));
  mod.fill_th(C.th, W.data(), p[6], p[7]);
  S = C;
  double cur_ll  = mod.loglik(C);
  double cur_lpr = logprior(p.data(), has_xi, s2_mu, s2_lam, s2_xi,
                            a_alpha, b_alpha, a_sig, b_sig);

  // Blocks, with the caches each one invalidates:
  //   0 = (muq, lamq, lalpha)   -> gamma margin
  //   1 = (muy, lamy, lsigma)   -> normal margin
  //   2 = (xi0, xi1)            -> copula parameter
  //   3 = (lalpha, xi0, xi1)    -> gamma margin and copula parameter
  //
  // Block 3 is the remedy for the dominant mixing problem. The shape alpha and
  // the copula parameters are strongly dependent in the posterior, because
  // residual dependence between the two markers can be absorbed either by a
  // tight gamma margin (large alpha) or by a strong copula; the fitted values
  // run from alpha = 66 under Gumbel to alpha = 269 under independence. Blocks
  // 0 and 2 can only move along the coordinate axes of that ridge, so neither
  // travels along it. Appending a joint update of (lalpha, xi0, xi1), whose
  // adapted covariance learns the ridge, lets the chain move along it directly.
  // It costs one extra gamma evaluation per iteration.
  const int bidx[4][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, -1}, {2, 6, 7}};
  const int blen[4] = {3, 3, 2, 3};
  const bool inv_q[4]  = {true,  false, false, true};
  const bool inv_y[4]  = {false, true,  false, false};
  const bool inv_th[4] = {false, false, true,  true};
  const int nblock = has_xi ? (cross_block ? 4 : 3) : 2;
  std::vector<std::vector<double> > lam(4), mu(4), acc(4);
  std::vector<std::vector<double> > Sig(4);
  for (int b = 0; b < 4; ++b) {
    lam[b].assign(blen[b], 2.38 * 2.38 / blen[b]);
    acc[b].assign(blen[b], 0.0);
    mu[b].assign(blen[b], 0.0);
    Sig[b].assign(blen[b] * blen[b], 0.0);
    for (int j = 0; j < blen[b]; ++j) {
      mu[b][j] = p[bidx[b][j]];
      Sig[b][j * blen[b] + j] = 0.01;
    }
  }
  std::vector<double> lamW(n, 1.0), muW(W), varW(n, 1.0), accW(n, 0.0);

  NumericMatrix PAR(8 + pd, M);
  std::vector<double> Wsum(n, 0.0);
  int Wcnt = 0;
  const int nWs = (thin_W > 0) ? (M - burn) / thin_W : 0;
  NumericMatrix Wsamp(nWs > 0 ? n : 0, nWs > 0 ? nWs : 0);
  int wj = 0;

  // online criteria states (Welford), for joint / q|y / y|q
  std::vector<double> mp[3], ml[3], vl[3], mip[3];
  for (int k = 0; k < 3; ++k) {
    mp[k].assign(n, 0.0); ml[k].assign(n, 0.0);
    vl[k].assign(n, 0.0); mip[k].assign(n, 0.0);
  }
  int ccnt = 0;

  std::vector<double> lpW(n), propW(n), lppW(n), aWv(n);
  RNGScope scope;

  for (int t = 1; t <= M; ++t) {
    const double gam = 1.0 / std::sqrt((double)t);

    // ---- Steps 1-3: parameter blocks ------------------------------------
    for (int b = 0; b < nblock; ++b) {
      const int L = blen[b];
      std::vector<double> sdb(L);
      for (int j = 0; j < L; ++j) {
        const double v = Sig[b][j * L + j];
        sdb[j] = std::sqrt(v > 1e-12 ? v : 1e-12);
      }
      // A joint proposal costs one gamma evaluation per block instead of L, but
      // it is a different algorithm: it must use the correlations in the adapted
      // covariance (via its Cholesky factor, ridged for positive definiteness)
      // and the L-dimensional optimal acceptance rate rather than the scalar one.
      std::vector<double> ch;
      if (blockwise && L > 1) {
        ch.assign(L * L, 0.0);
        for (int r = 0; r < L; ++r) {
          for (int c2 = 0; c2 <= r; ++c2) {
            double s = Sig[b][r * L + c2];
            if (r == c2) s += 1e-10;
            for (int k = 0; k < c2; ++k) s -= ch[r * L + k] * ch[c2 * L + k];
            if (r == c2) ch[r * L + c2] = std::sqrt(s > 1e-12 ? s : 1e-12);
            else         ch[r * L + c2] = s / ch[c2 * L + c2];
          }
        }
      }
      const double tgt = (blockwise && L > 1) ? BLOCK_TARGET[L] : TARGET;

      const int npro = blockwise ? 1 : L;
      for (int j = 0; j < npro; ++j) {
        std::vector<double> pp(p);
        if (blockwise) {
          std::vector<double> zz(L);
          for (int k = 0; k < L; ++k) zz[k] = norm_rand();
          for (int r = 0; r < L; ++r) {
            double e = 0.0;
            if (L > 1) { for (int c2 = 0; c2 <= r; ++c2) e += ch[r * L + c2] * zz[c2]; }
            else e = sdb[0] * zz[0];
            pp[bidx[b][r]] = p[bidx[b][r]] + std::sqrt(lam[b][0]) * e;
          }
        } else {
          pp[bidx[b][j]] = p[bidx[b][j]] + std::sqrt(lam[b][j]) * sdb[j] * norm_rand();
        }

        double a = 0.0;
        if (pp[4] >= 0.0) {               // identifiability: lam_y >= 0
          // recompute only the pieces this block touches
          if (inv_q[b])
            mod.fill_q(S.lf1, S.u1, S.lu1, S.z1, W.data(), pp[0], pp[1], std::exp(pp[2]));
          if (inv_y[b])
            mod.fill_y(S.lf2, S.u2, S.lu2, S.z2, W.data(), pp[3], pp[4], std::exp(pp[5]));
          if (inv_th[b])
            mod.fill_th(S.th, W.data(), pp[6], pp[7]);

          const double prop_ll  = mod.loglik(S);
          const double prop_lpr = logprior(pp.data(), has_xi, s2_mu, s2_lam, s2_xi,
                                           a_alpha, b_alpha, a_sig, b_sig);
          const double lr = (prop_ll + prop_lpr) - (cur_ll + cur_lpr);
          a = R_finite(lr) ? (lr >= 0.0 ? 1.0 : std::exp(lr)) : 0.0;
          if (unif_rand() < a) {
            p = pp; cur_ll = prop_ll; cur_lpr = prop_lpr;
            if (inv_q[b])  { C.lf1 = S.lf1; C.u1 = S.u1; C.lu1 = S.lu1; C.z1 = S.z1; }
            if (inv_y[b])  { C.lf2 = S.lf2; C.u2 = S.u2; C.lu2 = S.lu2; C.z2 = S.z2; }
            if (inv_th[b]) C.th = S.th;
            if (blockwise) { for (int k = 0; k < L; ++k) acc[b][k] += 1.0; }
            else acc[b][j] += 1.0;
          } else {
            if (inv_q[b])  { S.lf1 = C.lf1; S.u1 = C.u1; S.lu1 = C.lu1; S.z1 = C.z1; }
            if (inv_y[b])  { S.lf2 = C.lf2; S.u2 = C.u2; S.lu2 = C.lu2; S.z2 = C.z2; }
            if (inv_th[b]) S.th = C.th;
          }
        }
        if (blockwise) { for (int k = 0; k < L; ++k)
            lam[b][k] = std::exp(std::log(lam[b][k]) + gam * (a - tgt)); }
        else lam[b][j] = std::exp(std::log(lam[b][j]) + gam * (a - tgt));
      }

      std::vector<double> d(L);
      for (int k = 0; k < L; ++k) {
        d[k] = p[bidx[b][k]] - mu[b][k];
        mu[b][k] += gam * d[k];
      }
      for (int r = 0; r < L; ++r)
        for (int c2 = 0; c2 < L; ++c2)
          Sig[b][r * L + c2] += gam * (d[r] * d[c2] - Sig[b][r * L + c2]);
    }

    // ---- Step 4: beta, exact conjugate draw ------------------------------
    std::vector<double> XtW(pd, 0.0);
    for (int k = 0; k < pd; ++k) {
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += X(i, k) * W[i];
      XtW[k] = s;
    }
    std::vector<double> zb(pd);
    for (int k = 0; k < pd; ++k) zb[k] = norm_rand();
    for (int k = 0; k < pd; ++k) {
      double m = 0.0;
      for (int l = 0; l < pd; ++l) m += betapost_cov(k, l) * XtW[l];
      double e = 0.0;
      for (int l = 0; l < pd; ++l) e += betapost_cholT(k, l) * zb[l];
      beta[k] = m + e;
    }

    // ---- Step 5: latent variables ----------------------------------------
    std::vector<double> Xb(n, 0.0);
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int k = 0; k < pd; ++k) s += X(i, k) * beta[k];
      Xb[i] = s;
    }
    for (int i = 0; i < n; ++i) {
      const double d = W[i] - Xb[i];
      lpW[i] = mod.ld(C, i) - 0.5 * d * d - LOG_SQRT_2PI;
    }
    for (int i = 0; i < n; ++i) {
      const double v = varW[i] > 1e-8 ? varW[i] : 1e-8;
      propW[i] = W[i] + std::sqrt(lamW[i] * v) * norm_rand();
    }
    mod.fill_q(S.lf1, S.u1, S.lu1, S.z1, propW.data(), p[0], p[1], std::exp(p[2]));
    mod.fill_y(S.lf2, S.u2, S.lu2, S.z2, propW.data(), p[3], p[4], std::exp(p[5]));
    mod.fill_th(S.th, propW.data(), p[6], p[7]);
    for (int i = 0; i < n; ++i) {
      const double d = propW[i] - Xb[i];
      lppW[i] = mod.ld(S, i) - 0.5 * d * d - LOG_SQRT_2PI;
      const double lr = lppW[i] - lpW[i];
      double a = R_finite(lr) ? (lr >= 0.0 ? 1.0 : std::exp(lr)) : 0.0;
      if (!R_finite(a)) a = 0.0;
      aWv[i] = a;
    }
    for (int i = 0; i < n; ++i) {
      if (unif_rand() < aWv[i]) {
        W[i] = propW[i];
        C.lf1[i] = S.lf1[i]; C.u1[i] = S.u1[i]; C.lu1[i] = S.lu1[i]; C.z1[i] = S.z1[i];
        C.lf2[i] = S.lf2[i]; C.u2[i] = S.u2[i]; C.lu2[i] = S.lu2[i]; C.z2[i] = S.z2[i];
        C.th[i]  = S.th[i];
        accW[i] += 1.0;
      } else {
        S.lf1[i] = C.lf1[i]; S.u1[i] = C.u1[i]; S.lu1[i] = C.lu1[i]; S.z1[i] = C.z1[i];
        S.lf2[i] = C.lf2[i]; S.u2[i] = C.u2[i]; S.lu2[i] = C.lu2[i]; S.z2[i] = C.z2[i];
        S.th[i]  = C.th[i];
      }
      lamW[i] = std::exp(std::log(lamW[i]) + gam * (aWv[i] - TARGET));
      const double mo = muW[i];
      muW[i]  = mo + gam * (W[i] - mo);
      varW[i] = varW[i] + gam * ((W[i] - mo) * (W[i] - mo) - varW[i]);
    }
    cur_ll = mod.loglik(C);

    for (int k = 0; k < 8; ++k)  PAR(k, t - 1) = p[k];
    for (int k = 0; k < pd; ++k) PAR(8 + k, t - 1) = beta[k];

    // ---- criteria ---------------------------------------------------------
    if (t > burn) {
      if (track_criteria) {
        ++ccnt;
        for (int i = 0; i < n; ++i) {
          const double lc = (family == INDEP) ? 0.0
            : log_cop(C.u1[i], C.u2[i], C.lu1[i], C.lu2[i], C.z1[i], C.z2[i],
                      C.th[i], family);
          const double v[3] = { C.lf1[i] + C.lf2[i] + lc,
                                C.lf1[i] + lc,
                                C.lf2[i] + lc };
          for (int k = 0; k < 3; ++k) {
            const double lv = R_finite(v[k]) ? v[k] : R_NegInf;
            mp[k][i]  += (std::exp(lv) - mp[k][i]) / ccnt;
            mip[k][i] += (std::exp(-lv) - mip[k][i]) / ccnt;
            const double d1 = lv - ml[k][i];
            ml[k][i] += d1 / ccnt;
            vl[k][i] += d1 * (lv - ml[k][i]);
          }
        }
      }
      for (int i = 0; i < n; ++i) Wsum[i] += W[i];
      ++Wcnt;
      if (thin_W > 0 && (t - burn) % thin_W == 0 && wj < nWs) {
        for (int i = 0; i < n; ++i) Wsamp(i, wj) = W[i];
        ++wj;
      }
    }

    if (verbose > 0 && t % verbose == 0)
      Rcpp::Rcout << "  [" << family << "] " << t << "/" << M
                  << "  loglik " << cur_ll << "\n";
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();
  }

  NumericVector Wmean(n);
  for (int i = 0; i < n; ++i) Wmean[i] = Wsum[i] / (Wcnt > 0 ? Wcnt : 1);

  NumericVector crit;
  if (track_criteria) {
    crit = NumericVector::create(
      _["WAIC"] = 0.0, _["CWAIC_qy"] = 0.0, _["CWAIC_yq"] = 0.0,
      _["CVML"] = 0.0, _["CCVML_qy"] = 0.0, _["CCVML_yq"] = 0.0);
    for (int k = 0; k < 3; ++k) {
      double lppd = 0.0, pen = 0.0, cv = 0.0;
      for (int i = 0; i < n; ++i) {
        lppd += std::log(mp[k][i] > 1e-300 ? mp[k][i] : 1e-300);
        pen  += vl[k][i] / ccnt;
        cv   -= std::log(mip[k][i] > 1e-300 ? mip[k][i] : 1e-300);
      }
      crit[k]     = -2.0 * lppd + 2.0 * pen;
      crit[3 + k] = cv;
    }
  }

  // acceptance is reported per block, since blocks 0 and 3 share lalpha
  NumericVector accb(nblock);
  CharacterVector accnm(nblock);
  const char* bnames[4] = {"marg_q", "marg_y", "copula", "alpha_xi"};
  for (int b = 0; b < nblock; ++b) {
    double s = 0.0;
    for (int j = 0; j < blen[b]; ++j) s += acc[b][j];
    accb[b] = s / (blen[b] * (double)M);
    accnm[b] = bnames[b];
  }
  accb.names() = accnm;
  double aw = 0.0;
  for (int i = 0; i < n; ++i) aw += accW[i];

  return List::create(
    _["PAR"] = PAR, _["W"] = NumericVector(W.begin(), W.end()),
    _["Wmean"] = Wmean, _["Wsamp"] = Wsamp,
    _["criteria"] = crit, _["acc_blocks"] = accb,
    _["acc_W"] = aw / (n * (double)M), _["M"] = M, _["burn"] = burn);
}
