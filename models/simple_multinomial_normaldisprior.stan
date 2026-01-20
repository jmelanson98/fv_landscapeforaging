// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=1> C;                // # colonies
  int<lower=1> K;                // # traps
  int<lower=1> O;                // # observations
  int starts[C];                 // start indices for each colony
  int lengths[C];                // lengths for each colony
  matrix[O, 2] trap_pos;         // trap coordinates
  vector[O] sample_effort;      // sampling effort for each obs
  int colony_id[O];              // colony id for each obs
  int trap_id[O];                // trap id for each obs
  int y_obs[O];                  // concatenated counts for all colonies
  vector[O] yn;                 // are there any bees from that colony at that trap? 0/1
  real lower_x;                 // hard bounds on colony locations
  real upper_x;                 // hard bounds on colony locations
  real lower_y;                 // hard bounds on colony locations
  real upper_y;                 // hard bounds on colony locations
  real Rmax;
}


parameters {
  real<lower=0> rho;
  real<lower=0> sigma;
  vector[K] eps;

  // using hard bounds here so that init doesn't fail with - inf
  vector<lower=lower_x, upper=upper_x>[C] delta_x;
  vector<lower=lower_y, upper=upper_y>[C] delta_y;
}

transformed parameters {
  // epsilon scaling
  vector[K] eps_scale = eps*sqrt(sigma);
}

model {
  // set priors
  rho ~ lognormal(log(0.5), 0.5);
  sigma ~ exponential(1);
  eps ~ normal(0, 1);
  
  {
  
  // calculate likelihood
  for (c in 1:C) {
    
    // get some indices
    int start = starts[c]; // index start
    int length = lengths[c]; // index length
    
    // put prior on distance to traps where we find bees
    vector[length] dis = sqrt( square(delta_x[colony_id[start]] - trap_pos[start:start+length-1,1]) +
                       square(delta_y[colony_id[start]] - trap_pos[start:start+length-1,2]) );
    for (l in 1:length){
      if(yn[start+l-1] == 1){
        target += normal_lpdf(dis[l] | 0, Rmax/2.58);
      }
    }
    
    // compute visitation intensity lambda for each trap in that landscape
    vector[length] lambda_ik = -0.5*(dis / rho)^2 + 
                      eps_scale[trap_id[start:start+length-1]] + 
                      log(sample_effort[start:start+length-1]);
    
    // add to target likelihood
    y_obs[start:start+length-1] ~ multinomial(softmax(lambda_ik));
  }
}
}

