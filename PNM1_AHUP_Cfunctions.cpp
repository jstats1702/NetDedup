#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double logp_Xi (const double &I, const uvec &Xi) 
{
     // Hierarchical Uniform Prior (HUP)
     double r2 = I - max(Xi), r1 = I - 2 * r2, n;
     double out = 0.0;
     if (r2 == 0) {
          for (n = 2; n <= I; ++n) out += log(n);
     } else {
          if (r2 == floor(I/2)) {
               out += r2 * log(2.0);
               for (n = 2; n <= r2; ++n) out += log(n);
          } else {
               out += r2 * log(2.0); 
               for (n = 2; n <= r2; ++n) out += log(n);
               for (n = 2; n <= r1; ++n) out += log(n);
          }
     }
     return( out );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &pi_new, const rowvec &u_new, const uword &i, const double &I, const double &L, const double &K, const uvec &Mell, const vec &Vartheta, const umat &W, uvec Xi, mat Pi, const mat &Pe, const double &sigsq, const double &beta, mat Ustar, const uvec &Y)
{
     // computes the log conditional
     Xi(i) = Xi_i;
     relabel(Xi);
     if (caso == 1 && Xi_i != curlab) {
          Pi.shed_row(curlab-1);
          Ustar.shed_row(curlab-1);
     }
     if (Xi_i == newlab) {
          Pi = join_vert(Pi, pi_new);
          Ustar = join_vert(Ustar, u_new);
     }
     double out = ll_iter_pro_i(i, I, L, Mell, Vartheta, W, Xi, Pi, Pe) + ll_iter_net_i(i, I, beta, Ustar, Xi, Y);
     for (uword l = 0; l < L; ++l) 
          out += log(Vartheta( accu(Mell.rows(0, l)) - Mell(l) + Pi(Xi(i)-1, l) - 1 ));
     for (uword k = 0; k < K; ++k) 
          out += R::dnorm(Ustar(Xi(i)-1, k), 0.0, sqrt(sigsq), 1);
     out += logp_Xi(I, Xi);
     return( out );
}

void sample_Xi (const double &I, const double &L, const double &K, const uvec &Mell, const vec &Psi, const vec &Vartheta, umat &W, uvec &Xi, mat &Pi, const mat &Pe, const double &sigsq, const double &beta, mat &Ustar, const uvec &Y)
{
     // option 2 : Alg. 5 Neal(2000)
     uword i, l, curlab, newlab, prolab, caso, newval;
     double accp;
     rowvec pi_new(L), u_new(K);
     for (i = 0; i < I; ++i) {
          curlab = Xi(i);        // current label
          newlab = max(Xi) + 1;  // new label
          u_new  = sqrt(sigsq) * randn<rowvec>(K);
          for (l = 0; l < L; ++l) 
               pi_new(l) = wsample(Vartheta.rows( accu(Mell.rows(0, l)) - Mell(l), accu(Mell.rows(0, l)) - 1 )) + 1;
          // support newval
          uvec val = Xi;
          get_val_two(i, newlab, L, W, val, Pi, Pe);
          // draw proposal
          vec probs(val.n_elem, fill::ones); 
          prolab = val( wsample(probs) );  // proposal
          uvec id = find(Xi == curlab);    // id: actor indices pointing to curlab 
          caso = id.n_elem;                // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          // compute acceptance rate
          if (prolab == curlab)
               accp = 0.0;
          else
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, u_new, i, I, L, K, Mell, Vartheta, W, Xi, Pi, Pe, sigsq, beta, Ustar, Y) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, u_new, i, I, L, K, Mell, Vartheta, W, Xi, Pi, Pe, sigsq, beta, Ustar, Y) );
          // set new label
          if (R::runif(0.0, 1.0) < accp) 
               newval = prolab;
          else 
               newval = curlab;
          // modify Xi, Pi and W accordingly
          Xi(i) = newval;
          relabel(Xi);
          if (caso == 1 && newval != curlab) {
               Pi.shed_row(curlab-1);
               Ustar.shed_row(curlab-1);
          }
          if (newval == newlab) {
               Pi = join_vert(Pi, pi_new);
               sample_W_i(i, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
               Ustar = join_vert(Ustar, u_new);
          }
     }
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const uvec &Mell, const mat &Pe, const double &K, const uvec &Y, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path)
{
     double L = Mell.n_elem;  // n fields
     double I = Xi.n_elem;  // n subjects  
     // hyperparameters
     vec alpha(accu(Mell), fill::ones);
     double a_sig = 2.0 + pow(0.5, -2.0); // CV0 = 0.5
     double b_sig = (a_sig - 1.0) * ( sqrt(I)/(sqrt(I) - 2.0) ) * ( ( pow(datum::pi, K/2.0) / exp( lgamma(K/2.0 + 1.0) ) ) * pow(I, 2.0/K) );
     double omesq = 100.0;     
     // parameter initialization
     vec  Psi(L);
     vec  Vartheta(accu(Mell));
     umat W(I, L);
     mat  Pi(max(Xi), L);
     prior_int_pro(I, L, a_psi, b_psi, Mell, alpha, Psi, Vartheta, W, Xi, Pi, Pe);
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
          sample_Xi       (I, L, K, Mell, Psi, Vartheta, W, Xi, Pi, Pe, sigsq, beta, Ustar, Y);
          sample_Pi       (L, Mell, Vartheta, W, Xi, Pi, Pe);
          sample_W        (I, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_U        (s, nstar_N, nstar_U, delta_U, I, K, sigsq, beta, Ustar, Xi, Y);
          sample_Vartheta (L, Mell, alpha, Vartheta, W, Pi, Pe);
          sample_Psi      (I, L, a_psi, b_psi, Psi, W);
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
               Rcout << " *** PROFILE and NETWORK Model ***  " << 100.0 * myround(s/S, 3) << "% complete"
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