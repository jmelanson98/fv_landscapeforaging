//Pope & Jha, 2017, Cons. Genetics
// Modified to use rho rather than beta parametrization of length scale

data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
matrix[K, 2] trap; // trap coordinates 
int y[C, K]; // counts of bees in traps
real lowerbound; // uniform prior on colony location 
real upperbound; // uniform prior on colony location
vector[K] floral; // floral quality at traps
real<lower=0> priorVa; // prior variance on std deviations
real<lower=0> priorCo; // prior variance on other coefficients
real<lower=0> rho_center; //prior median for length parameter, !! on log scale!!
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

// calculate intensity
for(k in 1:K){
  for(i in 1:C){
    dis[i, k] = sqrt(square(delta[i, 1] - trap[k,1]) + square(delta[i, 2] - trap[k,2]));
    lambda[i, k] = dis[i,k]/(rho*exp(theta*floral[k])) + mu + zeta_scale[i] + eps_scale[k];
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
        lambda[i, k] = dis[i,k]/(rho*exp(theta*floral[k])) + mu + zeta_scale[i] + eps_scale[k];
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

