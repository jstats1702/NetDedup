#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double logp_Xi (const double &I, const double &a, const double &q, const double &eta, const double &theta_NB, const uvec &Xi) 
{
     double N = max(Xi);
     double out = lgamma(N + a) + N * ( log(q) + eta * log(1.0 - theta_NB) - log(1.0 - pow(1.0 - theta_NB, eta)) );
     for (uword n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          out += lgamma(id.n_elem + eta) - lgamma(eta);
     }
     return( out );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &u_new, const uword &i, const double &I, const double &K, const double &a, const double &q, const double &eta, const double &theta_NB, uvec Xi, const double &sigsq, const double &beta, mat Ustar, const uvec &Y)
{
     // computes the log conditional
     Xi(i) = Xi_i;
     relabel(Xi);
     if (caso == 1 && Xi_i != curlab) {
          Ustar.shed_row(curlab-1);
     }
     if (Xi_i == newlab) {
          Ustar = join_vert(Ustar, u_new);
     }
     double out = ll_iter_net_i(i, I, beta, Ustar, Xi, Y);
     for (uword k = 0; k < K; ++k) 
          out += R::dnorm(Ustar(Xi(i)-1, k), 0.0, sqrt(sigsq), 1);
     out += logp_Xi(I, a, q, eta, theta_NB, Xi);
     return( out );
}


void sample_Xi (const double &I, const double &K, const double &a, const double &q, const double &eta, const double &theta_NB, uvec &Xi, const double &sigsq, const double &beta, mat &Ustar, const uvec &Y)
{
     // option 2 : Alg. 5 Neal(2000)
     uword i, curlab, newlab, prolab, caso, newval;
     double accp;
     rowvec u_new(K);
     for (i = 0; i < I; ++i) {
          curlab = Xi(i);        // current label
          newlab = max(Xi) + 1;  // new label
          u_new  = sqrt(sigsq) * randn<rowvec>(K);
          // support newval
          uvec val = Xi;
          get_val_any_net(i, newlab, val);
          // draw proposal
          vec probs(val.n_elem, fill::ones); 
          prolab = val( wsample(probs) );  // proposal
          uvec id = find(Xi == curlab);    // id: actor indices pointing to curlab 
          caso = id.n_elem;                // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          // compute acceptance rate
          if (prolab == curlab)
               accp = 0.0;
          else
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, u_new, i, I, K, a, q, eta, theta_NB, Xi, sigsq, beta, Ustar, Y) - logcond_Xi(curlab, caso, curlab, newlab, u_new, i, I, K, a, q, eta, theta_NB, Xi, sigsq, beta, Ustar, Y) );
          // set new label
          if (R::runif(0.0, 1.0) < accp) 
               newval = prolab;
          else 
               newval = curlab;
          // modify Xi and Ustar accordingly
          Xi(i) = newval;
          relabel(Xi);
          if (caso == 1 && newval != curlab) {
               Ustar.shed_row(curlab-1);
          }
          if (newval == newlab) {
               Ustar = join_vert(Ustar, u_new);
          }
     }
}

double logcond_eta (const double &I, const double &a_eta, const double &b_eta, const double &eta, const double &theta_NB, const uvec &Xi)
{
     double N = max(Xi);
     double out = (a_eta - 1.0) * log(eta) - b_eta * eta + N * ( eta * log(1.0 - theta_NB) - log(1.0 - pow(1.0 - theta_NB, eta)) );
     for (uword n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          out += lgamma(id.n_elem + eta) - lgamma(eta);
     }
     return( out );
}

void sample_eta (const double &s, double &mean_c_eta, double &var_c_eta, double &nstar_eta, double &delta_eta, const double &I, const double &a_eta, const double &b_eta, double &eta, const double &theta_NB, const uvec &Xi)
{
     // eta = r
     // Metropolis step
     double x = log(eta);
     double x_p;
     if (R::runif(0.0, 1.0) < 0.95)
          x_p = R::rnorm(x, delta_eta * sqrt(var_c_eta));
     else 
          x_p = R::rnorm(x, 0.1);
     double r = exp( logcond_eta(I, a_eta, b_eta, exp(x_p), theta_NB, Xi) + x_p - logcond_eta(I, a_eta, b_eta, exp(x), theta_NB, Xi) - x );
     if (R::runif(0.0, 1.0) < r) {
          x = x_p; 
          ++nstar_eta;
     }
     eta = exp(x);
     tunning(delta_eta, 0.44, 0.025, nstar_eta/s, s);
     update_moments_scarlar(s, mean_c_eta, var_c_eta, eta);
}

double logcond_theta_NB (const double &I, const double &a_t, const double &b_t, const double &eta, const double &theta_NB, const uvec &Xi)
{
     double N = max(Xi);
     return ( (a_t + I - 1.0) * log(theta_NB) + (b_t + eta * N - 1.0) * log(1.0 - theta_NB) - N * log(1.0 - pow(1.0 - theta_NB, eta)) );
}

void sample_theta_NB (const double &s, double &mean_c_theta_NB, double &var_c_theta_NB, double &nstar_theta_NB, double &delta_theta_NB, const double &I, const double &a_t, const double &b_t, const double &eta, double &theta_NB, const uvec &Xi)
{
     // theta_NB = p
     // Metropolis step
     double x = log(theta_NB) - log(1.0 - theta_NB);
     double x_p;
     if (R::runif(0.0, 1.0) < 0.95)
          x_p = R::rnorm(x, delta_theta_NB * sqrt(var_c_theta_NB));
     else 
          x_p = R::rnorm(x, 0.1);
     double r = exp( logcond_theta_NB(I, a_t, b_t, eta, 1.0/(1.0 + exp(-x_p)), Xi) + ( x_p - 2.0 * log(1.0 + exp(x_p)) ) - logcond_theta_NB(I, a_t, b_t, eta, 1.0/(1.0 + exp(-x)), Xi) - ( x - 2.0 * log(1.0 + exp(x)) ) );
     if (R::runif(0.0, 1.0) < r) {
          x = x_p; 
          ++nstar_theta_NB;
     }
     theta_NB = 1.0/(1.0 + exp(-x));
     tunning(delta_theta_NB, 0.44, 0.025, nstar_theta_NB/s, s);
     update_moments_scarlar(s, mean_c_theta_NB, var_c_theta_NB, theta_NB);
}

//[[Rcpp::export]]
void MCMC (uvec Xi, const double &a, const double &q, const double &a_eta, const double &b_eta, const double &a_t, const double &b_t, const double &K, const uvec &Y, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path)
{
     double I = Xi.n_elem;  // n subjects  
     // hyperparameters
     double a_sig = 2.0 + pow(0.5, -2.0); // CV0 = 0.5
     double b_sig = (a_sig - 1.0) * ( sqrt(I)/(sqrt(I) - 2.0) ) * ( ( pow(datum::pi, K/2.0) / exp( lgamma(K/2.0 + 1.0) ) ) * pow(I, 2.0/K) );
     double omesq = 100.0;     
     // parameter initialization
     double eta = R::rgamma(a_eta, 1.0/b_eta);
     double theta_NB = R::rbeta(a_t, b_t);
     double sigsq;
     double beta;
     mat Ustar(max(Xi), K);
     prior_init_net(a_sig, b_sig, sigsq, omesq, beta, Ustar);
     // Metropolis
     double nstar_eta = 0.0;
     double delta_eta = 2.38/sqrt(1.0);
     double mean_c_eta = eta;
     double var_c_eta  = 0.0;
     double nstar_theta_NB = 0.0;
     double delta_theta_NB = 2.38/sqrt(1.0);
     double mean_c_theta_NB = theta_NB;
     double var_c_theta_NB  = 0.0;
     double nstar_N = 0.0;
     double nstar_U = 0.0;
     double delta_U = 2.38/sqrt(K);
     double nstar_beta = 0.0;
     double delta_beta = 2.38/sqrt(1.0);
     double mean_c_beta = beta;
     double var_c_beta  = 0.0;
     // write samples: opening files
     char* full;
     string nam;
     nam = "ll_chain";    full = mypaste0(path, nam); ofstream ll_chain;    ll_chain.open(full);
     nam = "Xi_chain";    full = mypaste0(path, nam); ofstream Xi_chain;    Xi_chain.open(full);
     nam = "DIC1";        full = mypaste0(path, nam); ofstream DIC1;        DIC1.open(full);
     nam = "DIC2";        full = mypaste0(path, nam); ofstream DIC2;        DIC2.open(full);
     nam = "WAIC1";       full = mypaste0(path, nam); ofstream WAIC1;       WAIC1.open(full);
     nam = "WAIC2";       full = mypaste0(path, nam); ofstream WAIC2;       WAIC2.open(full);
     nam = "Theta_hat";   full = mypaste0(path, nam); ofstream Theta_hat;   Theta_hat.open(full);
     nam = "sigsq_chain"; full = mypaste0(path, nam); ofstream sigsq_chain; sigsq_chain.open(full);
     nam = "beta_chain";  full = mypaste0(path, nam); ofstream beta_chain;  beta_chain.open(full);
     // information criteria parameters
     double ll;
     double ll_mean = 0.0;
     double theta;
     vec lppd(Y.n_elem, fill::zeros);
     vec mell(Y.n_elem, fill::zeros);
     vec sqll(Y.n_elem, fill::zeros);
     vec That(Y.n_elem, fill::zeros);
     // chain
     double S = nburn + nskip * nsams, s;
     uword i, ii, idx;
     for (s = 1; s <= S; ++s) {
          // update
          sample_Xi       (I, K, a, q, eta, theta_NB, Xi, sigsq, beta, Ustar, Y);
          sample_U        (s, nstar_N, nstar_U, delta_U, I, K, sigsq, beta, Ustar, Xi, Y);
          sample_eta      (s, mean_c_eta, var_c_eta, nstar_eta, delta_eta, I, a_eta, b_eta, eta, theta_NB, Xi);
          sample_theta_NB (s, mean_c_theta_NB, var_c_theta_NB, nstar_theta_NB, delta_theta_NB, I, a_t, b_t, eta, theta_NB, Xi);
          sample_beta     (s, mean_c_beta, var_c_beta, nstar_beta, delta_beta, I, omesq, beta, Ustar, Xi, Y);
          sample_sigsq    (K, a_sig, b_sig, sigsq, Ustar, Xi);
          // save samples
          if ((uword)s > nburn && (uword)s % nskip == 0) {
               for (i = 0; i < I; ++i) Xi_chain << Xi(i) << " ";
               Xi_chain << "\n";
               ll = ll_iter_net(I, beta, Ustar, Xi, Y);
               ll_chain    << ll    << "\n";
               sigsq_chain << sigsq << "\n";
               beta_chain  << beta  << "\n";
               // compute information criteria inputs
               ll_mean += ll/S;
               for (i = 1; i < I; ++i) {
                    for (ii = 0; ii < i; ++ii) {
                         idx = dyx(i, ii, I);
                         theta = 1.0/(1.0 + exp( -( beta - norm( Ustar.row( Xi(i) - 1 ) - Ustar.row( Xi(ii) - 1 ) ) ) ) );
                         That(idx) += theta/S;
                         if (Y(idx) == 1) {
                              lppd(idx) += theta/S;
                              mell(idx) += log(theta)/S;
                              sqll(idx) += pow(log(theta), 2.0);
                         } else {
                              lppd(idx) += (1.0 - theta)/S;
                              mell(idx) += (1.0 - log(theta))/S;
                              sqll(idx) += pow(log(1.0 - theta), 2.0);
                         }
                    }
               }
          }
          // display
          if ((uword)s % ndisp == 0) {
               vec r = myallelic(Xi);
               Rcout << " *** NETWORK Model ***  " << 100.0 * myround(s/S, 3) << "% complete"
                     << ", mr_eta = "   << myround(nstar_eta/s, 3)
                     << ", mr_theta = " << myround(nstar_theta_NB/s, 3)
                     << ", mr_beta = "  << myround(nstar_beta/s, 3)
                     << ", mr_U = "     << myround(mean(nstar_U/nstar_N), 3)
                     << ", N = "        << max(Xi)
                     << ", SR = "       << myround(r(0)/I, 4) 
                     << endl;
          }
     }  // end chain
     // compute information criteria
     double ll_hat = 0.0;  // log p(y|theta hat)
     for (i = 1; i < I; ++i) {
          for (ii = 0; ii < i; ++ii) {
               idx = dyx(i, ii, I);
               if (Y(idx) == 1)
                    ll_hat += log(That(idx));
               else
                    ll_hat += log(1.0 - That(idx));
          }
     }
     // save information criteria and Theta hat
     DIC1  <<  2.0 * ll_hat            - 2.0 * ( 2.0 * (ll_hat - ll_mean) )                                  << "\n";
     WAIC1 <<  2.0 * accu(log( lppd )) - 2.0 * ( 2.0 * ( accu(log( lppd )) - ll_mean ) )                     << "\n";
     DIC2  << -2.0 * ll_hat            + 2.0 * ( 2.0 * ( accu( ( sqll - S * pow(mell, 2.0) )/(S - 1.0) ) ) ) << "\n";
     WAIC2 << -2.0 * accu(log( lppd )) + 2.0 * ( 1.0 * ( accu( ( sqll - S * pow(mell, 2.0) )/(S - 1.0) ) ) ) << "\n";  // yes! it should be 1.0
     for (i = 0; i < That.n_elem; ++i) Theta_hat << That(i) << "\n";
     // close files
     ll_chain.close();
     Xi_chain.close();
     DIC1.close();
     DIC2.close();
     WAIC1.close();
     WAIC2.close();
     Theta_hat.close();
     sigsq_chain.close();
     beta_chain.close();
}