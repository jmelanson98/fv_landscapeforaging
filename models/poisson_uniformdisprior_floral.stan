// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=1> C;                // # colonies
  int<lower=1> K;                // # traps
  int<lower=1> O;                // # observations
  matrix[O, 2] trap_pos;         // trap coordinates
  vector[O] sample_effort;      // sampling effort for each obs
  int colony_id[O];              // colony id for each obs
  int trap_id[O];                // trap id for each obs
  vector[O] fq;                  // floral resources at the trap
  int y_obs[O];                  // concatenated counts for all colonies
  vector[O] yn;                 // are there any bees from that colony at that trap? 0/1
  real lower_x;                 // hard bounds on colony locations
  real upper_x;                 // hard bounds on colony locations
  real lower_y;                 // hard bounds on colony locations
  real upper_y;                 // hard bounds on colony locations
  real Rmax;
  real steepness;
}


parameters {
  real<lower=0> rho;
  real<lower=0> sigma;
  real theta;
  real alpha;
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
  theta ~ normal(0,1);
  alpha ~ normal(0,3);
  sigma ~ exponential(1);
  eps ~ normal(0, 1);
  
  {
  
  // calculate likelihood
  for (i in 1:O) {
    
    // put prior on distance to traps where we find bees
    real dis = sqrt( square(delta_x[colony_id[i]] - trap_pos[i,1]) +
                       square(delta_y[colony_id[i]] - trap_pos[i,2]) );
    if(yn[i] == 1){
        target += -log1p_exp((dis-Rmax)*steepness);
      }
    
    
    // compute visitation intensity lambda for each trap in that landscape
    real lambda_ik = alpha-0.5*(dis / rho)^2 + theta*fq[i] +
                      eps_scale[trap_id[i]] + 
                      log(sample_effort[i]);
    
    // add to target likelihood
    y_obs[i] ~ poisson_log(lambda_ik);
  }
}
}
