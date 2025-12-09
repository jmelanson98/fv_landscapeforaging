// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=1> C;                // # colonies
  int<lower=1> K;                // # traps
  int<lower=1> O;                // # observations
  int starts[C];                 // start indices for each colony
  int lengths[C];                // lengths for each colony
  matrix[O, 2] trap_pos;         // trap coordinates
  vector[O] sample_effort;          // sampling effort for each obs
  int colony_id[O];              // colony id for each obs
  int trap_id[O];                // trap id for each obs
  int y_obs[O];                  // concatenated counts for all colonies
  vector[O] yn;                     // are there any bees from that colony at that trap? 0/1
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
  vector[K] eps;
  
  // using hard bounds here so that init doesn't fail with - inf
  array[C] real<lower=lower_x, upper=upper_x> delta_x;
  array[C] real<lower=lower_y, upper=upper_y> delta_y;
}

model {
  // set priors
  rho ~ lognormal(log(0.5), 0.5);
  alpha ~ normal(0,1);
  sigma ~ normal(0, 1);
  eps ~ normal(0, 1);
  
  {
  // get eps_scale
  vector[K] eps_scale = eps*sqrt(sigma);
  
  // calculate likelihood
  for (c in 1:C) {
    
    int start = starts[c]; // index start
    int length = lengths[c]; // index length
    
    // get data subsets
    int y_ik[length] = y_obs[start : start+length-1]; // yobs for colony i
    vector[length] yn_ik = yn[start : start+length-1]; // was any bee observed for ik?
    vector[length] sample_effort_ik = sample_effort[start:start+length-1];
    int k[length] = trap_id[start:start+length-1]; // trap indices
    int i = colony_id[start]; // colony id
    matrix[length,2] trap_i = trap_pos[start:start+length-1,];
    
    // compute lambda for each trap in that landscape
    vector[length] dis = sqrt( square(delta_x[i] - trap_i[,1]) +
                       square(delta_y[i] - trap_i[,2]) );
    vector[length] lambda_ik = alpha - 0.5*(dis / rho)^2 + eps_scale[k] + log(sample_effort_ik);
    
    // compute multinomial probabilities and add to target likelihood
    vector[length] multi_probs = softmax(lambda_ik);
    y_ik ~ multinomial(multi_probs);

    // penalize distances that are outside Rmax
    // a bit like a soft indicator function, to maintain differentiability for HMC
    target += yn_ik .* (-penalty * log1p_exp((dis-Rmax)*steepness));
      
    }
  }
}
