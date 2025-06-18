// Parallelized exponential model using reduce_sum to speed up computation over many colonies
// J. Melanson (with help from chatgpt)
// June 17, 2025

functions {
  real partial_log_lik(array[] int slice_idx,
                       int start, int end,
                       int K,
                       matrix trap,
                       array[,] int y,
                       vector floral,
                       real rho,
                       real sigma_sqrt,
                       real tau_sqrt,
                       real mu,
                       real theta,
                       vector eps,
                       vector zeta,
                       array[] delta) {
                       
    real total = 0;
    for (i in start:end) {
      for (k in 1:K) {
        real dis = sqrt(square(delta[i, 1] - trap[k,1]) + square(delta[i, 2] - trap[k,2]));
        real lambda = dis/(-rho * exp(theta * floral[k])) + mu + zeta[i]*tau_sqrt + eps[k]*sigma_sqrt;
        total += poisson_log_lpmf(y[i,k] | lambda);
      }
    }
    return total;
  }
}


data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
matrix[K, 2] trap; // trap coordinates 
array[C,K] int y; // counts of bees in traps
real lowerbound; // uniform prior on colony location 
real upperbound; // uniform prior on colony location
vector[K] floral; // floral quality at traps
real<lower=0> priorVa; // prior variance on std deviations
real<lower=0> priorCo; // prior variance on other coefficients
real<lower=0> rho_center; //prior median for length parameter, !! on log scale!!
real<lower=0> rho_sd; //prior sd of length parameter, !!on log scale!!
int<lower=1> grainsize; // number of colonies to parallelize per chunk
}

transformed data {
int colony_ids[C];
for (i in 1:C)
  colony_ids[i] = i;
}

parameters { // see text for definitions
real<lower=0> rho; 
real<lower=0> sigma;
real<lower=0> tau; 
real theta;
real mu; 
vector[K] eps; 
vector[C] zeta;
array [C,2] vector<lower=lowerbound, upper=upperbound> delta;

}


transformed parameters {
  real<lower=0> tau_sqrt = sqrt(tau);
  real<lower=0> sigma_sqrt = sqrt(sigma); 
  vector[C] zeta_scale = zeta * tau_sqrt;
  vector[K] eps_scale = eps*sigma_sqrt;
}
  
  
model {
  
// temporary declarations
matrix[C,K] dis; //distance from colony C to trap K?
matrix[C,K] lambda; //rate of captures for colony C at trap K?

// priors
sigma ~ normal(0, priorVa);
tau ~ normal(0, priorVa); 
rho ~ lognormal(rho_center, rho_sd);
mu ~ normal(0, priorCo); 
theta ~ normal(0, priorCo);

// random effects for traps
eps ~ normal(0, 1); 
zeta ~ normal(0, 1);

// calculate intensity using partial_log_lik function
target += reduce_sum(partial_log_lik, colony_ids, grainsize,
                     K, trap, y, floral,
                     rho, sigma_sqrt, tau_sqrt, mu, theta,
                     eps, zeta, delta);

}

// remove GQ block for now -- can't use reduce_sum there
// do these calculations outside of Stan, OR reintroduce GQ block using map_rect (very confusing...)
