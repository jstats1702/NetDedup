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

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &pi_new, const uword &i, const double &I, const double &L, const uvec &Mell, const double &vt, const vec &Vartheta, umat W, uvec Xi, mat Pi, const mat &Pe)
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
     out += logp_Xi(I, vt, Xi);
     return( out );
}

void sample_Xi (const double &I, const double &L, const uvec &Mell, const double &vt, const vec &Psi, const vec &Vartheta, umat &W, uvec &Xi, mat &Pi, const mat &Pe)
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
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, i, I, L, Mell, vt, Vartheta, W, Xi, Pi, Pe) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, i, I, L, Mell, vt, Vartheta, W, Xi, Pi, Pe) );
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

void sample_vt (const double &I, const double &a_vt, const double &b_vt, double &vt, const uvec &Xi)
{
     double r2 = I - max(Xi);
     vt = R::rbeta(a_vt + r2, b_vt + floor(I/2.0) - r2);
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const double &a_vt, const double &b_vt, const uvec &Mell, const mat &Pe, const uword &nburn, const uword &nsams, const uword &nskip, const uword& ndisp, string path)
{
     double L = Mell.n_elem;  // n fields
     double I = Pe.n_rows;    // n subjects
     // hyperparameters
     vec alpha(accu(Mell), fill::ones);
     // parameter initialization
     double vt = R::rbeta(a_vt, b_vt);
     vec  Psi(L);
     vec  Vartheta(accu(Mell));
     umat W(I, L);
     mat  Pi(max(Xi), L);
     prior_int_pro(I, L, a_psi, b_psi, Mell, alpha, Psi, Vartheta, W, Xi, Pi, Pe);
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
          sample_Xi       (I, L, Mell, vt, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_Pi       (L, Mell, Vartheta, W, Xi, Pi, Pe);
          sample_W        (I, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_vt       (I, a_vt, b_vt, vt, Xi);
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