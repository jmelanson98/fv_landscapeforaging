//STAN code from Pope & Jha, 2017, Cons. Genetics

data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
matrix[K, 2] trap; // trap coordinates 
int y[C, K]; // counts of bees in traps
real lowerbound; // uniform prior on colony location 
real upperbound; // uniform prior on colony location
vector[K] floral; // floral quality at traps
real<lower=0> priorVa; // prior variance on std deviations 
real<lower=0> priorBe; // prior variance on beta
real<lower=0> priorCo; // prior variance on other coefficients
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
  eps_scale = eps*sigma_sqrt; //multiply eps by standard deviation of dist it was drawn from?? (because eps prior is N(0,1) and the true variance is not 1, it's eps)...is this to make the model run better? like a more efficient parametrization?
  zeta_scale = zeta*tau_sqrt; //multiply zeta by standard deviation of dist it was drawn from??
  
  }
  
  
model {
  
// temporary declarations
matrix[C,K] dis; //distance from colony C to trap K?
matrix[C,K] lambda; //rate of captures for colony C at trap K?

// priors
sigma ~ normal(0, priorVa);
tau ~ normal(0, priorVa); 
beta ~ normal(0, priorBe);
mu ~ normal(0, priorCo); 
theta ~ normal(0, priorCo);

// random effects for traps
eps ~ normal(0, 1); 
zeta ~ normal(0, 1);

// calculate intensity
for(k in 1:K){
  for(i in 1:C){
    dis[i, k] = sqrt(square(delta[i][1] - trap[k,1]) + square(delta[i][2] - trap[k,2]));
    lambda[i, k] = -beta*dis[i, k] + theta*floral[k] + mu + zeta_scale[i] + eps_scale[k];
    y[i, k] ~ poisson_log(lambda[i, k]);
  } 
}

}  
  
// log likelihood
// this bit changed...
//increment_log_prob( -exp(lambda) ); //increment_log_prob: deprecated? use: target += exp(lambda); ?
//increment_log_prob( y .* lambda ); }
// because i'm not that fancy
// also added the likelihood calculation to the above for-loop, because Stan does not like taking a matrix in poisson_log()
