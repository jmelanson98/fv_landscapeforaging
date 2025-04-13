//Pope & Jha, 2017, Cons. Genetics
// Modified to include ypred in generated quantities block

data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
matrix[K, 2] trap; // trap coordinates 
int y[C, K]; // counts of bees in traps
real lowerbound; // uniform prior on colony location 
real upperbound; // uniform prior on colony location
vector[K] floral; // floral quality at traps
real<lower=0> priorVa; // prior variance on coefficients 
real<lower=0> priorCo; // prior variance on std deviations
}

parameters { // see text for definitions
real<lower=0> beta; 
real<lower=0> sigma;
real<lower=0> tau; 
real theta;
real mu; 
vector[K] eps; 
vector[C] zeta;
// the following is no longer valid stan syntax:
// matrix <lower=lowerbound , upper=upperbound >[C,2] delta; }
// change to:
array [C] vector<lower=lowerbound, upper=upperbound>[2] delta;

}


transformed parameters { 
  real<lower=0> tau_sqrt;
  real<lower=0> sigma_sqrt; 
  vector[C] zeta_scale;
  vector[K] eps_scale;
  
  tau_sqrt = sqrt(tau); 
  sigma_sqrt = sqrt(sigma); 
  eps_scale = eps*sigma_sqrt; //non-centered parametrization
  zeta_scale = zeta*tau_sqrt; //non-centered parametrization
  
  // temporary declarations
matrix[C,K] dis; //distance from colony C to trap K?
matrix[C,K] lambda; //rate of captures for colony C at trap K?

  // distance and lambda
for(k in 1:K){
  for(i in 1:C){
    dis[i, k] = sqrt(square(delta[i][1] - trap[k,1]) + square(delta[i][2] - trap[k,2]));
    lambda[i, k] = -beta*dis[i, k] + theta*floral[k] + mu + zeta_scale[i] + eps_scale[k];
  } 
}
}
  
  
model {

// priors
sigma ~ normal(0, priorVa);
tau ~ normal(0, priorVa); 
beta ~ normal(0, priorCo);
mu ~ normal(0, priorCo); 
theta ~ normal(0, priorCo);

// random effects for traps
eps ~ normal(0, 1); 
zeta ~ normal(0, 1);

// calculate intensity
for(k in 1:K){
  for(i in 1:C){
    y[i, k] ~ poisson_log(lambda[i, k]);
  } 
}
}

// will this work? who knows!?
generated quantities {
  int y_rep[C, K];        // Posterior predictive samples
  for (k in 1:K) {
    for (i in 1:C){
      y_rep[i, k] = poisson_log_rng(lambda[i, k]);  // Simulate new data points
    }
  }
}

