// Estimate foraging distance from Bombus sib captures
// May 23, 2025
// J Melanson

// Modified from Pope & Jha, 2017, Cons. Genetics
// New features:
// exponentiated quadratic distance decay
// informative prior on colony location
functions {
  real bilinear_interp(matrix landcover_mat, real x, real y) {
    int x0 = to_int(floor(x));
    int x1 = x0 + 1;
    int y0 = to_int(floor(y));
    int y1 = y0 + 1;

    // check that we're within bounds of matrix
    if (x0 < 1 || x1 > cols(landcover_mat) || y0 < 1 || y1 > rows(landcover_mat))
      return negative_infinity(); // outside matrix

    real dx = x - x0;
    real dy = y - y0;

    // get values at the 4 cells surrounding proposed colony location
    real Q11 = landcover_mat[y0, x0];
    real Q21 = landcover_mat[y0, x1];
    real Q12 = landcover_mat[y1, x0];
    real Q22 = landcover_mat[y1, x1];

    // bilinear interpolation
    // i.e., interpolate first in the x direction, then y direction
    // (this function is what you get after a bit of algebraic fiddling)
    return (1 - dx) * (1 - dy) * Q11 +
           dx * (1 - dy) * Q21 +
           (1 - dx) * dy * Q12 +
           dx * dy * Q22;
  }
}


data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
int<lower=1> L; // landscape size
int<lower=1> res; //resolution of nesting landscape matrix

matrix[K, 2] trap; // trap coordinates
vector[K] fq; // floral quality at traps
int y[C, K]; // counts of bees in traps

matrix[L, L] nesting_landscape; // matrix of nesting quality values for study landscape
real lowerbound; // lower limit on colony location 
real upperbound; // upper limit on colony location

real<lower=0> priorVa; // prior variance on std deviations
real<lower=0> priorCo; // prior variance on other coefficients
real<lower=0> rho_center; //prior median for length parameter !!on log scale!!
real<lower=0> rho_sd; //prior sd of length parameter, !!on log scale!!
}

parameters { // see text for definitions
real<lower=0> rho; 
real<lower=0> sigma;
real<lower=0> tau; 
real theta;
real mu; 
vector[K] eps; 
vector[C] zeta;
array [C] vector<lower=lowerbound, upper=upperbound>[2] delta;

}


transformed parameters {
  real<lower=0> tau_sqrt = sqrt(tau);
  real<lower=0> sigma_sqrt = sqrt(sigma); 
  vector[C] zeta_scale = zeta * tau_sqrt;
  vector[K] eps_scale = eps*sigma_sqrt;
}
  
  
model {
  
// temporary declarations
matrix[C,K] dis; //distance from colony C to trap K
matrix[C,K] lambda; //rate of captures for colony C at trap K

// priors
sigma ~ normal(0, priorVa);
tau ~ normal(0, priorVa); 
rho ~ lognormal(rho_center, rho_sd);  // log-normal with median = rho_center
mu ~ normal(0, priorCo); 
theta ~ normal(0, priorCo);

// random effects for traps
eps ~ normal(0, 1); 
zeta ~ normal(0, 1);

for(i in 1:C){
  // custom prior for colony location
  real suitability = bilinear_interp(nesting_landscape, delta[i,1], delta[i,2]);
  target += suitability;
  
  for(k in 1:K){
    // condition on observed data
    dis[i, k] = sqrt(square(delta[i, 1] - trap[k,1]) + square(delta[i, 2] - trap[k,2]));
    lambda[i, k] = -0.5*(dis[i,k]/rho)^2 + theta*fq[k] + mu + zeta_scale[i] + eps_scale[k];
    y[i, k] ~ poisson_log(lambda[i, k]);
  } 
}
}

generated quantities {
  vector[C] colony_dist;        // Declare estimated colony foraging distance
  
  colony_dist = rep_vector(0, C);  // Initialize colony distance to 0
  
  // start anonymous scope
  {
    matrix[C,K] dis; //distance from colony C to trap K
    matrix[C,K] lambda; //rate of captures for colony C at trap K
    vector[C] V;     // Declare local variable for each colony's total visitation
    real alpha = 1e-12;  // set smaaaaall number to prevent division by 0
    
    // Recompute lambda and dis
    for(k in 1:K){
      for(i in 1:C){
        dis[i, k] = sqrt(square(delta[i, 1] - trap[k,1]) + square(delta[i, 2] - trap[k,2]));
        lambda[i, k] = -0.5*(dis[i,k]/rho)^2 + theta*fq[k] + mu + zeta_scale[i] + eps_scale[k];
      } 
    }
    
    // Compute V for normalization
    for (i in 1:C){
      V[i] = sum(exp(lambda[i,]));
    }
    
    
    // compute colony_dist to be saved outside the anonymous scope
    for (k in 1:K){
      colony_dist = colony_dist + (dis[,k] .* exp(lambda[,k]) ./ (V + alpha));
    }
  }
}

