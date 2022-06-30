#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

vec myallelic3 (const uvec &x)
{
     uword M = 3;
     vec nk = mytable(x), out(M);
     for (uword i = 0; i < M; ++i) {
          uvec id = find(nk == (i+1));
          out(i) = id.n_elem;
     }
     return( out );
}

double logfactorial (const double &x)
{
     double out = 0.0;
     if ( x > 1.0 ) {
          for (uword n = 2; n <= x; ++n) out += log((double)n);     
     }
     return( out );
}

void nothree (uvec &x) 
{
     // removes pairs from x
     uvec ux = unique(x);
     for (uword i = 0; i < ux.n_elem; ++i) {
          uvec rep = find(x == ux(i));
          if (rep.n_elem > 2) {
               for (uword k = 0; k < rep.n_elem; ++k) x.shed_row(rep(k) - k);
          }
     }
     x = unique(x); 
}

void get_val_three (const uword &i, const uword &newlab, const double &L, const umat &W, uvec &val, const mat &Pi, const mat &Pe)
{
     val.shed_row(i);
     nothree(val);  // remove all "pairs" from val
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

double logp_Xi (const double &I, const double &vt2, const double &vt3, const uvec &Xi) 
{
     vec r = myallelic3(Xi);
     double M  = 3.0;
     double Q3 = floor(I/M);
     double Q2 = floor((I - 3*r(2))/2.0);
     return( r(1) * log(2.0) + r(2) * log(6.0) + logfactorial(r(0)) + logfactorial(r(1)) + logfactorial(r(2)) + R::dbinom(r(2), Q3, vt3, 1) + R::dbinom(r(1), Q2, vt2, 1) );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &pi_new, const uword &i, const double &I, const double &L, const uvec &Mell, const double &vt2, const double &vt3, const vec &Vartheta, umat W, uvec Xi, mat Pi, const mat &Pe)
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
     out += logp_Xi(I, vt2, vt3, Xi);
     return( out );
}

void sample_Xi (const double &I, const double &L, const uvec &Mell, const double &vt2, const double &vt3, const vec &Psi, const vec &Vartheta, umat &W, uvec &Xi, mat &Pi, const mat &Pe)
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
          get_val_three(i, newlab, L, W, val, Pi, Pe);
          // draw proposal
          vec probs(val.n_elem, fill::ones); 
          prolab = val( wsample(probs) );  // proposal
          uvec id = find(Xi == curlab);    // id: actor indices pointing to curlab 
          caso = id.n_elem;                // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          // compute acceptance rate
          if (prolab == curlab)
               accp = 0.0;
          else
               accp = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, i, I, L, Mell, vt2, vt3, Vartheta, W, Xi, Pi, Pe) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, i, I, L, Mell, vt2, vt3, Vartheta, W, Xi, Pi, Pe) );
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

void sample_vt (const double &I, const double &a_vt2, const double &b_vt2, const double &a_vt3, const double &b_vt3, double &vt2, double &vt3, const uvec &Xi)
{
     vec r = myallelic3(Xi);
     double M  = 3.0;
     double Q3 = floor(I/M);
     double Q2 = floor((I - 3*r(2))/2.0);
     vt2 = R::rbeta(a_vt2 + r(1), b_vt2 + Q2 - r(1));
     vt3 = R::rbeta(a_vt3 + r(2), b_vt3 + Q3 - r(2));
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const double &a_vt2, const double &b_vt2, const double &a_vt3, const double &b_vt3, const uvec &Mell, const mat &Pe, const uword &nburn, const uword &nsams, const uword &nskip, const uword& ndisp, string path)
{
     double L = Mell.n_elem;  // n fields
     double I = Pe.n_rows;    // n subjects
     // hyperparameters
     vec alpha(accu(Mell), fill::ones);
     // parameter initialization
     double vt2 = R::rbeta(a_vt2, b_vt2);
     double vt3 = R::rbeta(a_vt3, b_vt3);
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
          sample_Xi       (I, L, Mell, vt2, vt3, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_Pi       (L, Mell, Vartheta, W, Xi, Pi, Pe);
          sample_W        (I, L, Mell, Psi, Vartheta, W, Xi, Pi, Pe);
          sample_vt       (I, a_vt2, b_vt2, a_vt3, b_vt3, vt2, vt3, Xi);
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