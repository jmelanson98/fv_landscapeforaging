// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=0> C;                // total number of colonies
  int<lower=0> K;                // total number of traps
  int<lower=0> O;                // total number of observations (observation for a colony at a trap, in a timepoint)
  vector[O] sample_effort;
  matrix[O,2] trap_pos;              // trap coordinates
  int colony_id[O];              // colony id for each observation
  int trap_id[O];                 // trap id for each observation
  int y_obs[O];                   // number of individuals observed
  real upper_y; // hard bound on colony locations
  real upper_x; // hard bound on colony locations
  real lower_y; // hard bound on colony locations
  real lower_x; // hard bound on colony locations
}

transformed data {
  real Rmax = 2.0;     // max distance of colonies from traps where they observed
  real steepness = 10;    // steepness of penalty function
  real penalty = 10;     // penalty scale
}

parameters {
  real<lower=0> rho;
  real alpha;
  real<lower=0> sigma;
  real<lower=0> phi;
  real<lower=0, upper=1> theta; // probability of drawing a zero (zero-inflation)
  vector[K] eps;
  
  // using hard bounds here so that init doesn't fail with - inf
  array[C] real<lower=lower_x, upper=upper_x> delta_x;
  array[C] real<lower=lower_y, upper=upper_y> delta_y;
}


transformed parameters {
  real<lower=0> sigma_sqrt = sqrt(sigma);
  vector[K] eps_scale = eps*sigma_sqrt;
}



model {
  // set priors
  rho ~ lognormal(log(0.5), 0.5);
  alpha ~ normal(0,1);
  sigma ~ normal(0, 1);
  eps ~ normal(0, 1);
  phi ~ lognormal(log(10), 1);
  // implicit uniform prior on theta
  
  // compute log(theta), log(1-theta) once per iteration
  real log_theta = log(theta);
  real log_not_theta = log(1-theta);
  
  // calculate likelihood
  for (n in 1:O){
   
    // compute eta
    real dis = sqrt( square(delta_x[colony_id[n]] - trap_pos[n,1]) +
                       square(delta_y[colony_id[n]] - trap_pos[n,2]) );
    real eta = alpha - 0.5*(dis / rho)^2 + eps_scale[trap_id[n]] + log(sample_effort[n]);
    
    // compute zinb probabilities and add to target likelihood
    if (y_obs[n] == 0) {
      target += log_sum_exp(
        log_theta,
        log_not_theta + neg_binomial_2_log_lpmf(0 | eta, phi)
      );
    } else {
      target += log_not_theta + neg_binomial_2_log_lpmf(y_obs[n] | eta, phi);
      
      // penalize distances outside Rmax
      target += -penalty * log1p_exp((dis-Rmax)*steepness);
    }
  }
}
