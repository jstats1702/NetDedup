#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double logp_Xi (const double &I, const double &a, const double &q, const vec &Mu, const uvec &Xi) 
{
     double N = max(Xi), n;
     double out = lgamma(N + a) + N * log(q);
     for (n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          out += lgamma(id.n_elem + 1.0) + log(Mu(id.n_elem - 1));
     }
     return( out );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &pi_new, const rowvec &u_new, const uword &i, const double &I, const double &L, const double &K, const uvec &Mell, const double &a, const double &q, const vec &Mu, const vec &Vartheta, const umat &W, uvec Xi, mat Pi, const mat &Pe, const double &sigsq, const double &beta, mat Ustar, const uvec &Y)
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
     out += logp_Xi(I, a, q, Mu, Xi);
     return( out );
}

void sample_Xi (const double &I, const double &L, const double &K, const uvec &Mell, const double &a, const double &q, const vec &Mu, const vec &Psi, const vec &Vartheta, umat &W, uvec &Xi, mat &Pi, const mat &Pe, const double &sigsq, const double &beta, mat &Ustar, const uvec &Y)
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
          get_val_any(i, newlab, L, W, val, Pi, Pe);
          // draw proposal
          vec probs(val.n_elem, fill::ones); 
          prolab = val( wsample(probs) );  // proposal
          uvec id = find(Xi == curlab);    // id: actor indices pointing to curlab 
          caso = id.n_elem;                // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          // compute acceptance rate
          if (prolab == curlab)
               accp = 0.0;
          else
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, u_new, i, I, L, K, Mell, a, q, Mu, Vartheta, W, Xi, Pi, Pe, sigsq, beta, Ustar, Y) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, u_new, i, I, L, K, Mell, a, q, Mu, Vartheta, W, Xi, Pi, Pe, sigsq, beta, Ustar, Y) );
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

void sample_Mu (const double &I, const double &al, const vec &Mu0, vec &Mu, const uvec &Xi)
{
     vec r = myallelic(Xi);
     r.insert_rows(I, 1);
     r(I) = 0.0;
     Mu = rdirichlet(al * Mu0 + r);
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const double &a, const double &q, const double &al, const vec &Mu0, const uvec &Mell, const mat &Pe, const double &K, const uvec &Y, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path)
{
     double L = Mell.n_elem;  // n fields
     double I = Xi.n_elem;    // n subjects  
     // hyperparameters
     vec alpha(accu(Mell), fill::ones);
     double a_sig = 2.0 + pow(0.5, -2.0); // CV0 = 0.5
     double b_sig = (a_sig - 1.0) * ( sqrt(I)/(sqrt(I) - 2.0) ) * ( ( pow(datum::pi, K/2.0) / exp( lgamma(K/2.0 + 1.0) ) ) * pow(I, 2.0/K) );
     double omesq = 100.0;     
     // parameter initialization
     vec  Mu = rdirichlet(al * Mu0);
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
     double nstar_N     = 0.0;
     double nstar_U     = 0.0;
     double delta_U     = 2.38/sqrt(K);
     double nstar_beta  = 0.0;
     double delta_beta  = 2.38/sqrt(1.0);
     double mean_c_beta = beta;
     double var_c_beta  = 0.0;
     // write samples: opening files
     char* full;
     string nam;
     nam = "inmat_out"; full = mypaste0(path, nam); ofstream inmat_out; inmat_out.open(full);
     nam = "stats_out"; full = mypaste0(path, nam); ofstream stats_out; stats_out.open(full);
     // outputs
     rowvec inmat(I*(I-1)/2, fill::zeros);
     // chain
     double S = nburn + nskip * nsams, s;
     for (s = 1; s <= S; ++s) {
          // update
          sample_Xi       (I, L, K, Mell, a, q, Mu, Psi, Vartheta, W, Xi, Pi, Pe, sigsq, beta, Ustar, Y);
          sample_Pi       (L, Mell, Vartheta, W, Xi, Pi, Pe);
          sample_W        (I, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_U        (s, nstar_N, nstar_U, delta_U, I, K, sigsq, beta, Ustar, Xi, Y);
          sample_Mu       (I, al, Mu0, Mu, Xi);
          sample_Vartheta (L, Mell, alpha, Vartheta, W, Pi, Pe);
          sample_Psi      (I, L, a_psi, b_psi, Psi, W);
          sample_beta     (s, mean_c_beta, var_c_beta, nstar_beta, delta_beta, I, omesq, beta, Ustar, Xi, Y);
          sample_sigsq    (K, a_sig, b_sig, sigsq, Ustar, Xi);
          // save samples
          if ((uword)s > nburn && (uword)s % nskip == 0) {
               inmat += get_inmat(Xi)/((double)nsams);
               vec r  = myallelic(Xi);
               stats_out << max(Xi) << " " << r(0) << " " << max(find(r)) + 1 << "\n";
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
     // save outputs
     for (uword i = 0; i < inmat.n_elem; ++i) inmat_out << inmat(i) << "\n";
     // close files
     stats_out.close();
     inmat_out.close();
}