#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double logp_Xi (const double &I, const double &vt, const uvec &Xi) 
{
     double r2 = I - max(Xi), r1 = I - 2 * r2, Q2 = floor(I/2.0), n;
     double out = 0.0;
     if (r2 == 0) {
          for (n = 2; n <= I; ++n) out += log(n);
     } else {
          if (r2 == Q2) {
               out += r2 * log(2.0);
               for (n = 2; n <= r2; ++n) out += log(n);
          } else {
               out += r2 * log(2.0); 
               for (n = 2; n <= r2; ++n) out += log(n);
               for (n = 2; n <= r1; ++n) out += log(n);
          }
     }
     out += R::dbinom(r2, Q2, vt, 1);
     return( out );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &u_new, const uword &i, const double &I, const double &K, const double &vt, uvec Xi, const double &sigsq, const double &beta, mat Ustar, const uvec &Y)
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
     out += logp_Xi(I, vt, Xi);
     return( out );
}

void sample_Xi (const double &I, const double &K, const double &vt, uvec &Xi, const double &sigsq, const double &beta, mat &Ustar, const uvec &Y)
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
          get_val_two_net(i, newlab, val);
          // draw proposal
          vec probs(val.n_elem, fill::ones); 
          prolab = val( wsample(probs) );  // proposal
          uvec id = find(Xi == curlab);    // id: actor indices pointing to curlab 
          caso = id.n_elem;                // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          // compute acceptance rate
          if (prolab == curlab)
               accp = 0.0;
          else
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, u_new, i, I, K, vt, Xi, sigsq, beta, Ustar, Y) - logcond_Xi(curlab, caso, curlab, newlab, u_new, i, I, K, vt, Xi, sigsq, beta, Ustar, Y) );
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

void sample_vt (const double &I, const double &a_vt, const double &b_vt, double &vt, const uvec &Xi)
{
     double r2 = I - max(Xi);
     vt = R::rbeta(a_vt + r2, b_vt + floor(I/2.0) - r2);
}

//[[Rcpp::export]]
void MCMC (uvec Xi, const double &a_vt, const double &b_vt, const double &K, const uvec &Y, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path)
{
     double I = Xi.n_elem;  // n subjects  
     // hyperparameters
     double a_sig = 2.0 + pow(0.5, -2.0); // CV0 = 0.5
     double b_sig = (a_sig - 1.0) * ( sqrt(I)/(sqrt(I) - 2.0) ) * ( ( pow(datum::pi, K/2.0) / exp( lgamma(K/2.0 + 1.0) ) ) * pow(I, 2.0/K) );
     double omesq = 100.0;     
     // parameter initialization
     double vt = R::rbeta(a_vt, b_vt);
     double sigsq;
     double beta;
     mat Ustar(max(Xi), K);
     prior_init_net(a_sig, b_sig, sigsq, omesq, beta, Ustar);
     // Metropolis
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
          sample_Xi    (I, K, vt, Xi, sigsq, beta, Ustar, Y);
          sample_U     (s, nstar_N, nstar_U, delta_U, I, K, sigsq, beta, Ustar, Xi, Y);
          sample_vt    (I, a_vt, b_vt, vt, Xi);
          sample_beta  (s, mean_c_beta, var_c_beta, nstar_beta, delta_beta, I, omesq, beta, Ustar, Xi, Y);
          sample_sigsq (K, a_sig, b_sig, sigsq, Ustar, Xi);
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