#ifndef samplers_H
#define samplers_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

vec rdirichlet (const vec &alpha) 
{
     uword K = alpha.n_elem;
     vec out(K);
     for (uword k = 0; k < K; ++k) out(k) = R::rgamma(alpha(k), 1.0);  // draw Gamma variables
     return( out/accu(out) );
}

uword wsample (vec probs) 
{
     // weighted sampling according to probs
     // return indices fron 0 to n - 1, where n = probs.n_elem
     probs = probs/accu(probs);
     double u = R::runif(0.0, 1.0);
     if (u < probs(0)) {
          return 0;
     } else {
          uword n = probs.n_elem, i, out = 0;
          vec probsum = cumsum(probs);
          for (i = 1; i < n; ++i) {  
               if ( (probsum(i-1) <= u) && (u < probsum(i)) ) {
                    out = i;
                    goto endloop;
               }
          }
          endloop:
               return( out );
     }
}

double myround (const double &value, const uword &digits)
{
     // otherwise it will return 'nan' due to the log10() of zero
     if (value == 0.0) {
          return 0.0;
     } else {
          double factor = pow(10.0, digits - ceil(log10(fabs(value))));
          return( round(value * factor)/factor );
     }
}

char* mypaste0 (string path, string name)
{
     stringstream strname;
     strname << path << name << ".txt";
     string fullname = strname.str();
     string::iterator p = fullname.begin();
     char* chr = &(*p);
     return( chr );
}

vec mytable (const uvec &x)
{
     // x has continuous labels from 1 to N 
     uword N = max(x);
     vec out(N);
     for (uword n = 0; n < N; ++n) {
          uvec id = find(x == (n+1));
          out(n) = id.n_elem;
     }
     return( out );
}

vec myallelic (const uvec &x)
{
     uword I = x.n_rows;
     vec nk = mytable(x), out(I);
     for (uword i = 0; i < I; ++i) {
          uvec id = find(nk == (i+1));
          out(i) = id.n_elem;
     }
     return( out );
}

uword dyx (const uword &i, const uword &ii, const double &I)
{
     // indices used to store values in and call values from Y and Z
     // (vectorized) linear index from a lower triangular matrix position
     return( I*(I - 1)/2 - (I - ii)*(I - ii - 1)/2 + i - ii - 1 );
}

double ll_iter_pro (const double &I, const double &L, const uvec &Mell, const vec &Vartheta, const umat &W, const uvec &Xi, const mat &Pi, const mat &Pe)
{
     // computes the log-likelihood in a given iteration of the MCMC
     uword i, l, lo_vt;
     double out = 0.0;
     for (l = 0; l < L; ++l) {
          lo_vt = accu(Mell.rows(0, l)) - Mell(l);
          for (i = 0; i < I; ++i) {
               if (W(i, l) == 1) {
                    out += log(Vartheta(lo_vt + Pe(i, l) - 1));
               } else {
                    if (Pe(i, l) != Pi(Xi(i) - 1, l)) {
                         out = -datum::inf;
                         goto endloop;
                    }
               }
          }
     }
     endloop:
          return( out );
}

double ll_iter_pro_i (const uword &i, const double &I, const double &L, const uvec &Mell, const vec &Vartheta, const umat &W, const uvec &Xi, const mat &Pi, const mat &Pe)
{
     // computes the log-likelihood in a given iteration of the MCMC
     uword l, lo_vt;
     double out = 0.0;
     for (l = 0; l < L; ++l) {
          lo_vt = accu(Mell.rows(0, l)) - Mell(l);
          if (W(i, l) == 1) {
               out += log(Vartheta(lo_vt + Pe(i, l) - 1));
          } else {
               if (Pe(i, l) != Pi(Xi(i) - 1, l)) {
                    out = -datum::inf;
                    goto endloop;
               }
          }
     }
     endloop:
          return( out );
}

double ll_iter_net (const double &I, const double &beta, const mat &Ustar, const uvec &Xi, const uvec &Y)
{
     // computes the log-likelihood in a given iteration of the MCMC
     uword i, ii;
     double eta, out = 0.0;
     for (i = 1; i < I; ++i) {  // lower triangualar matrix
          for (ii = 0; ii < i; ++ii) {
               eta = beta - norm( Ustar.row( Xi(i) - 1 ) - Ustar.row( Xi(ii) - 1 ) );
               out += eta * Y( dyx(i, ii, I) ) - log(1.0 + exp(eta));
          }
     }
     return( 2.0 * out );
}

double ll_iter_net_i (const uword &i, const double &I, const double &beta, const mat &Ustar, const uvec &Xi, const uvec &Y)
{
     // computes the contribution of the log-likelihood for a given i and j
     uword ii;
     double eta, out = 0.0;
     for (ii = 0; ii < i; ++ii) {
          eta = beta - norm( Ustar.row( Xi(i) - 1 ) - Ustar.row( Xi(ii) - 1 ) );
          out += eta * Y( dyx(i, ii, I) ) - log(1.0 + exp(eta));
     }
     for (ii = i+1; ii < I; ++ii) {
          eta = beta - norm( Ustar.row( Xi(ii) - 1 ) - Ustar.row( Xi(i) - 1 ) );
          out += eta * Y( dyx(ii, i, I ) ) - log(1.0 + exp(eta));
     }
     return( 2.0 * out );
}

void tunning (double &delta, const double &mr0, const double &epsilon, const double &mix_rate, const double &s)
{
     // tunning paramter calibration
     // ntun = 100
     if ( ((uword)s % 100 == 0) && (abs(mix_rate - mr0) > epsilon) ) {
          double tmp = delta, cont = 1.0;
          do {
               tmp = delta + (0.1/cont) * (mix_rate - mr0);
               ++cont;
          } while ( (tmp <= 0.0) && (cont <= 100.0) );
          if (tmp > 0) delta = tmp;
     }
}

void update_moments_scarlar (const double &n_c, double &mean_c, double &var_c, const double &obs_new) 
{
     // var computed dividing by n instead of n-1
     var_c  = (n_c * (var_c + pow(mean_c, 2.0)) + pow(obs_new, 2.0))/(n_c + 1.0) - pow((n_c * mean_c + obs_new)/(n_c + 1.0), 2.0);
     mean_c = (n_c * mean_c + obs_new)/(n_c + 1.0);
}

void sample_Pi (const double &L, const uvec &Mell, const vec &Vartheta, const umat &W, const uvec &Xi, mat &Pi, const mat &Pe)
{
     uword N = max(Xi), l, n;
     uvec lth_col(1);
     for (l = 0; l < L; ++l) {
          lth_col(0) = l;
          for (n = 0; n < N; ++n) {
               uvec Rnj = find(Xi == (n+1));
               uvec idx = find(W.submat(Rnj, lth_col) == 0);
               if (idx.n_elem > 0)
                    Pi(n, l) = Pe(Rnj(idx(0)), l);
               else
                    Pi(n, l) = wsample(Vartheta.rows( accu(Mell.rows(0, l)) - Mell(l), accu(Mell.rows(0, l)) - 1 )) + 1;
          }
     }
}

double logcond_U (const rowvec &u, const uword &n, const uvec &id, const double &I, const double &K, const double &sigsq, const double &beta, mat Ustar, const uvec &Xi, const uvec &Y) 
{
     Ustar.row(n) = u;
     double out = 0.0;
     for (uword k = 0; k < id.n_elem; ++k) out += ll_iter_net_i(id(k), I, beta, Ustar, Xi, Y);
     for (uword k = 0; k < K; ++k) out += R::dnorm(Ustar(n, k), 0.0, sqrt(sigsq), 1);
     return( out );
}

void sample_U (const double &s, double &nstar_N, double &nstar_U, double &delta_U, const double &I, const double &K, const double &sigsq, const double &beta, mat &Ustar, const uvec &Xi, const uvec &Y)
{
     // Xi = [\xi_1, ... , \xi_I] is indexed from 1 to N where N = max{Xi}
     // U  = [u*_{\xi_1}, ... , u*_{\xi_I}]^T : I x K matrix
     // U* = [u*_1, ... , u*_N]^T             : N x K matrix
     // Metropolis step
     uword N = max(Xi), n;
     rowvec u_p(K);
	double r;
     for (n = 0; n < N; ++n) {
          uvec id = find( Xi == (n+1) );
          u_p = Ustar.row(n) + delta_U * randn<rowvec>(K);
          r = exp( logcond_U(u_p, n, id, I, K, sigsq, beta, Ustar, Xi, Y) - logcond_U(Ustar.row(n), n, id, I, K, sigsq, beta, Ustar, Xi, Y) );
          if (R::runif(0.0, 1.0) < r) {
               Ustar.row(n) = u_p;
               ++nstar_U;
          }
     }
     nstar_N += N;
     tunning(delta_U, 0.37, 0.025, nstar_U/nstar_N, s);
}

void sample_W_i (const uword &i, const double &L, const uvec &Mell, const vec &Psi, const vec &Vartheta, umat &W, const uvec &Xi, const mat &Pi, const mat &Pe)
{
     uword l, lo_vt;
     double logprob;
     for (l = 0; l < L; ++l) {
          lo_vt = accu(Mell.rows(0, l)) - Mell(l);
          W(i, l) = 1;
          if (Pe(i, l) == Pi(Xi(i) - 1, l)) {
               logprob = log(Psi(l)) + log(Vartheta( lo_vt + Pe(i, l) - 1));
               logprob -= log(exp(logprob) + 1.0 - Psi(l));
               W(i, l) = R::rbinom(1.0, exp(logprob));
          }
     }
}

void sample_W (const double &I, const double &L, const uvec &Mell, const vec &Psi, const vec &Vartheta, umat &W, const uvec &Xi, const mat &Pi, const mat &Pe)
{
     uword i, l, lo_vt;
     double logprob;
     for (l = 0; l < L; ++l) {
          lo_vt = accu(Mell.rows(0, l)) - Mell(l);
          for (i = 0; i < I; ++i) {
               W(i, l) = 1;
               if (Pe(i, l) == Pi(Xi(i)-1, l)) {
                    logprob = log(Psi(l)) + log(Vartheta(lo_vt + Pe(i, l) - 1));
                    logprob -= log(exp(logprob) + 1.0 - Psi(l));
                    W(i, l) = R::rbinom(1.0, exp(logprob));
               }
          }
     }
}

void sample_Vartheta (const double &L, const uvec &Mell, const vec &alpha, vec &Vartheta, const umat &W, const mat &Pi, const mat &Pe)
{
     uword m, l, lo_vt;
     uvec lth_col(1);
     for (l = 0; l < L; ++l) {
          lo_vt = accu(Mell.rows(0, l)) - Mell(l);
          lth_col(0) = l;
          vec alpha_new(Mell(l));
          for (m = 0; m < Mell(l); ++m) {
               uvec idx_Pi  = find(Pi.col(l) == (m+1));
               uvec idx_Pe  = find(Pe.col(l) == (m+1));
               alpha_new(m) = alpha(lo_vt + m) + idx_Pi.n_elem + accu(W.submat(idx_Pe, lth_col));
          }
          Vartheta.rows( lo_vt, lo_vt + Mell(l) - 1 ) = rdirichlet(alpha_new);
     }
}

void sample_Psi (const double &I, const double &L, const double &a, const double &b, vec &Psi, const umat &W)
{
     double Suma;
     for (uword l = 0; l < L; ++l) {
          Suma = accu(W.col(l));
          Psi(l) = R::rbeta(a + Suma, b + I - Suma);
     }
}

double logcond_beta (const double &I, const double &omesq, const double &beta, const mat &Ustar, const uvec &Xi, const uvec &Y)
{
     return( ll_iter_net(I, beta, Ustar, Xi, Y) + R::dnorm(beta, 0.0, sqrt(omesq), 1) );
}

void sample_beta (const double &s, double &mean_c_beta, double &var_c_beta, double &nstar_beta, double &delta_beta, const double &I, const double &omesq, double &beta, const mat &Ustar, const uvec &Xi, const uvec &Y)
{
     // Metropolis step
     double beta_p;
     if (R::runif(0.0, 1.0) < 0.95)
          beta_p = R::rnorm(beta, delta_beta * sqrt(var_c_beta));
     else 
          beta_p = R::rnorm(beta, 0.1);
     double r = exp( logcond_beta(I, omesq, beta_p, Ustar, Xi, Y) - logcond_beta(I, omesq, beta, Ustar, Xi, Y) );
     if (R::runif(0.0, 1.0) < r) {
          beta = beta_p; 
          ++nstar_beta;
     }
     tunning(delta_beta, 0.44, 0.025, nstar_beta/s, s);
     update_moments_scarlar(s, mean_c_beta, var_c_beta, beta);
}

void sample_sigsq (const double &K, const double &a_sig, const double &b_sig, double &sigsq, const mat &Ustar, const uvec &Xi)
{
     double N = max(Xi);
     sigsq = 1.0/R::rgamma( a_sig + 0.5 * K * N, 1.0/(b_sig + 0.5 * accu(pow(Ustar, 2.0))) ) ;
}

void prior_int_pro (const double &I, const double &L, const double &a_psi, const double &b_psi, const uvec &Mell, const vec &alpha, vec &Psi, vec &Vartheta, umat &W, const uvec &Xi, mat &Pi, const mat &Pe)
{
     // parameter initialization
     uword N = max(Xi), i, n, l, lo_vt, up_vt;
     for (l = 0; l < L; ++l) {
          Psi(l) = R::rbeta(a_psi, b_psi);
          lo_vt = accu(Mell.rows(0, l)) - Mell(l);
          up_vt = lo_vt + Mell(l) - 1;
          Vartheta.rows(lo_vt, up_vt) = rdirichlet(alpha.rows(lo_vt, up_vt));
          for (n = 0; n < N; ++n) 
               Pi(n, l) = wsample(Vartheta.rows(lo_vt, up_vt)) + 1;
     }
     for (i = 0; i < I; ++i)
          sample_W_i (i, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
     
}

void prior_init_net (const double &a_sig, const double &b_sig, double &sigsq, const double &omesq, double &beta, mat &Ustar)
{
     sigsq = 1.0/R::rgamma(a_sig, 1.0/b_sig);
     beta = R::rnorm(0.0, sqrt(omesq));
     for (uword i = 0; i < Ustar.n_rows; ++i) {
          for (uword k = 0; k < Ustar.n_cols; ++k) {
               Ustar(i, k) = R::rnorm(0.0, sqrt(sigsq));
          }
     }
}

mat get_Ustar (const mat &U, const uvec &Xi)
{
     // create U* (N x K), where N is the number of unique latent individuals
     uword N = max(Xi), n;
     mat Ustar(N, U.n_cols);
     for (n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          Ustar.row(n) = U.row(id(0));
     }
     return( Ustar );
}

mat get_U (const mat &Ustar, const uvec &Xi) 
{
     // creates U (I x K), where I is the number of actors
     // Xi given with the "right" labels
     uword N = max(Xi), n, i;
     mat U(Xi.n_elem, Ustar.n_cols);
     for (n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          for (i = 0; i < id.n_elem; ++i) U.row(id(i)) = Ustar.row(n);
     }
     return( U );
}

void relabel (uvec &Xi) 
{
     // re-label latent individuals from 1 to N
     // N is the number of unique elements in Xi
     uvec uXi = unique(Xi);  // unique elements of Xi, sorted in ascending order
     uword N = uXi.n_elem;
     if (N < uXi(N-1)) {
          for (uword n = 0; n < N; ++n) {
               uvec id = find(Xi == uXi(n));
               for (uword i = 0; i < id.n_elem; ++i) Xi(id(i)) = n+1;
          }
     }
}

void norep (uvec &x) 
{
     // removes pairs from x
     uword cont = x.n_elem, nr, i = 0;
     do {
          uvec rep = find(x == x(i));
          nr = rep.n_elem;
          if (nr > 1) {
               for (uword k = 0; k < nr; ++k) x.shed_row(rep(k) - k);
               --i;
          }
          cont -= nr;
          ++i;
     } while (cont > 0);
}

void get_val_two (const uword &i, const uword &newlab, const double &L, const umat &W, uvec &val, const mat &Pi, const mat &Pe)
{
     val.shed_row(i);
     norep(val);  // remove all "pairs" from val
     // if w_ijl = 0 and x_ijl != y_cl, then Xi_ij = c is impossible
     uword nval0 = val.n_elem, c = 0, l, cond;
     for (uword k = 0; k < nval0; ++k) {
          l = 0;
          cond = 0;
          while ( l < L && cond == 0 ) {
               if (( W(i, l) == 0 ) && ( Pe(i, l) != Pi(val(k - c)-1, l) )) {
                    val.shed_row(k - c);
                    ++c;
                    cond = 1;
               }
               ++l;
          }
     }
     // add new cluster
     val.insert_rows(val.n_elem, 1);
     val(val.n_elem - 1) = newlab;
}

void get_val_any (const uword &i, const uword &newlab, const double &L, const umat &W, uvec &val, const mat &Pi, const mat &Pe)
{
     val.shed_row(i);
     val = unique(val);
     // if w_ijl = 0 and x_ijl != y_cl, then Xi_ij = c is impossible
     uword nval0 = val.n_elem, c = 0, l, cond;
     for (uword k = 0; k < nval0; ++k) {
          l = 0; cond = 0;
          while ( l < L && cond == 0 ) {
               if (( W(i, l) == 0 ) && ( Pe(i, l) != Pi(val(k - c)-1, l) )) {
                    val.shed_row(k - c);
                    ++c;
                    cond = 1;
               }
               ++l;
          }
     }
     // add new cluster
     val.insert_rows(val.n_elem, 1);
     val(val.n_elem - 1) = newlab;
}

void get_val_two_net (const uword &i, const uword &newlab, uvec &val)
{
     val.shed_row(i);
     norep(val);  // remove all "pairs" from val
     // add new cluster
     val.insert_rows(val.n_elem, 1);
     val(val.n_elem - 1) = newlab;
}

void get_val_any_net (const uword &i, const uword &newlab, uvec &val)
{
     val.shed_row(i);
     val = unique(val);
     // add new cluster
     val.insert_rows(val.n_elem, 1);
     val(val.n_elem - 1) = newlab;
}

rowvec get_inmat (const uvec &Xi) 
{
     // vec lower triangular portion of inmat
     uword I = Xi.n_elem, i, ii;
     rowvec A(I*(I - 1)/2, fill::zeros);
     for (i = 1; i < I; ++i) {
          for (ii = 0; ii < i; ++ii) {
               if ( Xi(i) == Xi(ii) ) A( dyx(i, ii, I) ) = 1;
          }
     }
     return( A );
}

#endif