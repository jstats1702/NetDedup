#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

uword wsample (const vec &probs) 
{
     // weighted sampling according to probs
     // return indices fron 0 to n - 1, where n = probs.n_elem
     double u = R::runif(0.0, 1.0);
     if (u < probs(0)) {
          return 0;
     } else {
          uword n = probs.n_elem, out = 0;
          vec probsum = cumsum(probs);
          for (uword i = 1; i < n; ++i) {  
               if ( (probsum(i-1) <= u) && (u < probsum(i)) ) {
                    out = i;
                    goto endloop;
               }
          }
          endloop:
               return( out );
     }
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

//[[Rcpp::export]]
uvec gen_Xi_NBNBP_C (const double &a, const double &q, const double &r, const double &p, const uword &M, const uword &I) 
{
     // Sampler for NBNB prior
     // Produces one sample from NBNB
     // Needs M gibbs samples bc sampling is not exact
     double beta = exp( log(q) + r*log(1.0 - p) - log(1.0 - pow(1.0 - p, r)) );
     uword i, j;
     uvec cluster(I), clusterA(I); for (i = 0; i < I; ++i) cluster(i) = i+1;
     for (j = 0; j < M; ++j) {
          clusterA = cluster;
          for (i = 0; i < I; ++i) {
               clusterA.shed_row(i);
               relabel(clusterA);
               vec probs = log(mytable(clusterA) + r);
               probs.insert_rows(probs.n_elem, 1);
               probs(probs.n_elem - 1) = log(max(clusterA) + a) + log(beta) + log(r);
               probs = exp(probs - max(probs));
               probs = probs/sum(probs);
               clusterA.insert_rows(i, 1);
               clusterA(i) = wsample(probs) + 1;
               cluster = clusterA;
          }
     }
     return( cluster );
}

//[[Rcpp::export]]
vec rdirichlet (const vec &alpha) 
{
     uword K = alpha.n_elem;
     vec out(K);
     for (uword k = 0; k < K; ++k) out(k) = R::rgamma(alpha(k), 1.0);  // draw Gamma variables
     return( out/accu(out) );
}

//[[Rcpp::export]]
uvec gen_Xi_NBDP_C (const double &a, const double &q, const double &al, const vec &Mu0, const uword &M, const uword &I) 
{
     // Sampler for NBNB prior
     // Produces one sample from NBNB
     // Needs M gibbs samples bc sampling is not exact
     uword i, j, k;
     uvec cluster(I), clusterA(I); for (i = 0; i < I; ++i) cluster(i) = i+1;
     vec Mu(I+1);
     for (j = 0; j < M; ++j) {
          clusterA = cluster;
          vec alphaDM = al * Mu0.rows(0, I-1) +  myallelic(cluster);
          alphaDM.insert_rows(I, 1);
          alphaDM(I) = al * Mu0(I);
          Mu = rdirichlet(alphaDM);
          for (i = 0; i < I; ++i) {
               clusterA.shed_row(i);
               relabel(clusterA);
               vec nk = mytable(clusterA);
               vec probs(nk.n_elem);
               for (k = 0; k < nk.n_elem; ++k) probs(k) = log(nk(k) + 1.0) + log(Mu(nk(k))) - log(Mu(nk(k)-1));
               probs.insert_rows(probs.n_elem, 1);
               probs(probs.n_elem - 1) = log(max(clusterA) + a) + log(q) + log(Mu(0));
               probs = exp(probs - max(probs));
               probs = probs/sum(probs);
               clusterA.insert_rows(i, 1);
               clusterA(i) = wsample(probs) + 1;
               cluster = clusterA;
          }
     }
     return( cluster );
}