#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double logp_Xi (const double &theta, const uvec &Xi) 
{
     double N = max(Xi), out = N * log(theta);
     for (uword n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          out += lgamma(id.n_elem);
     }
     return( out );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &pi_new, const uword &i, const double &I, const double &L, const uvec &Mell, const double &theta, const vec &Vartheta, umat W, uvec Xi, mat Pi, const mat &Pe)
{
     // computes the log conditional
     Xi(i) = Xi_i;
     relabel(Xi);
     if (caso == 1 && Xi_i != curlab) {
          Pi.shed_row(curlab-1);
     }
     if (Xi_i == newlab) {
          Pi = join_vert(Pi, pi_new);
     }
     double out = ll_iter_pro_i(i, I, L, Mell, Vartheta, W, Xi, Pi, Pe);
     for (uword l = 0; l < L; ++l) 
          out += log(Vartheta( accu(Mell.rows(0, l)) - Mell(l) + Pi(Xi(i)-1, l) - 1 ));
     out += logp_Xi(theta, Xi);
     return( out );
}

void sample_Xi (const double &I, const double &L, const uvec &Mell, const double &theta, const vec &Psi, const vec &Vartheta, umat &W, uvec &Xi, mat &Pi, const mat &Pe)
{
     // option 2 : Alg. 5 Neal(2000)
     uword i, l, curlab, newlab, prolab, caso, newval;
     double accp;
     rowvec pi_new(L);
     for (i = 0; i < I; ++i) {
          curlab = Xi(i);        // current label
          newlab = max(Xi) + 1;  // new label
          for (l = 0; l < L; ++l) 
               pi_new(l) = wsample(Vartheta.rows( accu(Mell.rows(0, l)) - Mell(l), accu(Mell.rows(0, l)) - 1 )) + 1;
          // support newval
          uvec val = Xi;
          get_val_any(i, newlab, L, W, val, Pi, Pe);
          // draw proposal
          uvec id = find(Xi == curlab);    // id: actor indices pointing to curlab 
          caso = id.n_elem;                // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          vec probs(val.n_elem, fill::ones); 
          prolab = val( wsample(probs) );  // proposal
          // compute acceptance rate
          if (prolab == curlab)
               accp = 0.0;
          else
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, i, I, L, Mell, theta, Vartheta, W, Xi, Pi, Pe) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, i, I, L, Mell, theta, Vartheta, W, Xi, Pi, Pe) );
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
          }
          if (newval == newlab) {
               Pi = join_vert(Pi, pi_new);
               sample_W_i(i, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
          }
     }
}

double logcond_theta (const double &a_t, const double &b_t, const double &theta, const uvec &Xi)
{
     double N = max(Xi), out = N * log(theta);
     for (uword n = 0; n < N; ++n) {
          uvec id = find(Xi == (n+1));
          out += lgamma(id.n_elem);
     }
     out += R::dgamma(theta, a_t, 1.0/b_t, 1);
     return( out );
}

void sample_theta (const double &s, double &mean_c_theta, double &var_c_theta, double &nstar_theta, double &delta_theta, const double &I, const double &a_t, const double &b_t, double &theta, const uvec &Xi)
{
     // Metropolis step
     double x = log(theta);
     double x_p;
     if (R::runif(0.0, 1.0) < 0.95)
          x_p = R::rnorm(x, delta_theta * sqrt(var_c_theta));
     else 
          x_p = R::rnorm(x, 0.1);
     double r = exp( logcond_theta(a_t, b_t, exp(x_p), Xi) + x_p - logcond_theta(a_t, b_t, exp(x), Xi) - x );
     if (R::runif(0.0, 1.0) < r) {
          x = x_p; 
          ++nstar_theta;
     }
     theta = exp(x);
     tunning(delta_theta, 0.44, 0.025, nstar_theta/s, s);
     update_moments_scarlar(s, mean_c_theta, var_c_theta, theta);
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const double &a_t, const double &b_t, const uvec &Mell, const mat &Pe, const uword &nburn, const uword &nsams, const uword &nskip, const uword& ndisp, string path)
{
     double L = Mell.n_elem;  // n fields
     double I = Pe.n_rows;    // n subjects
     // hyperparameters
     vec alpha(accu(Mell), fill::ones);
     // parameter initialization
     double theta = R::rgamma(a_t, b_t);
     vec  Psi(L);
     vec  Vartheta(accu(Mell));
     umat W(I, L);
     mat  Pi(max(Xi), L);
     prior_int_pro(I, L, a_psi, b_psi, Mell, alpha, Psi, Vartheta, W, Xi, Pi, Pe);
     // Metropolis
     double nstar_theta  = 0.0;
     double delta_theta  = 2.38/sqrt(1.0);
     double mean_c_theta = theta;
     double var_c_theta  = 0.0;
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
          sample_Xi       (I, L, Mell, theta, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_Pi       (L, Mell, Vartheta, W, Xi, Pi, Pe);
          sample_W        (I, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_theta    (s, mean_c_theta, var_c_theta, nstar_theta, delta_theta, I, a_t, b_t, theta, Xi);
          sample_Vartheta (L, Mell, alpha, Vartheta, W, Pi, Pe);
          sample_Psi      (I, L, a_psi, b_psi, Psi, W);
          // save samples
          if ((uword)s > nburn && (uword)s % nskip == 0) {
               inmat += get_inmat(Xi)/((double)nsams);
               vec r  = myallelic(Xi);
               stats_out << max(Xi) << " " << r(0) << " " << max(find(r)) + 1 << "\n";
          }
          // display
          if ((uword)s % ndisp == 0) {
               vec r = myallelic(Xi);
               Rcout << " *** PROFILE Model ***  " << 100.0 * myround(s/S, 3) << "% complete"
                     << ", N = "  << max(Xi)
                     << ", SR = " << myround(r(0)/I, 4) 
                     << endl;
          }
     }  // end chain
     // save outputs
     for (uword i = 0; i < inmat.n_elem; ++i) inmat_out << inmat(i) << "\n";
     // close files
     stats_out.close();
     inmat_out.close();
}