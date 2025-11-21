// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=1> C;                // # colonies
  int<lower=1> L;                // # landscapes
  int<lower=1> total_traps;      // # traps
  int<lower=0> total_obs;        // total number of observations (observation for a colony at a trap, not just trap observations)
  matrix[total_traps, 2] trap_pos; // trap coordinates
  int sample_effort[total_traps]; // per trap sampling effort
  int traps_start[L];            // 1-based index for trap_pos
  int traps_n[L];                // # traps per landscape
  int colony_land[C];            // site id for each colony
  int y_start[C];                // start index in y_flat for each colony
  int y_n[C];                    // number of traps (length) for each colony's y
  int y_flat[total_obs];                  // concatenated counts for all colonies
  real upper_y; // hard bound on colony locations
  real upper_x; // hard bound on colony locations
  real lower_y; // hard bound on colony locations
  real lower_x; // hard bound on colony locations
}

transformed data {
  real Rmax = 2.0;     // max distance of colonies from traps where they observed
  real steepness = 50;    // steepness of penalty function
  real penalty = 2000;     // penalty scale
}

parameters {
  real<lower=0> rho; 
  real<lower=0> beta;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[total_traps] eps; 
  vector[C] zeta;
  
  // using hard bounds here so that init doesn't fail with - inf
  array[C] real<lower=lower_x, upper=upper_x> delta_x;
  array[C] real<lower=lower_y, upper=upper_y> delta_y;
}


transformed parameters {
  real<lower=0> tau_sqrt = sqrt(tau);
  real<lower=0> sigma_sqrt = sqrt(sigma); 
  vector[C] zeta_scale = zeta * tau_sqrt;
  vector[total_traps] eps_scale = eps*sigma_sqrt;
}

model {
  // set priors
  // see inside for loop for delta_x, delta_y priors
  rho ~ lognormal(log(0.5), 0.5);
  beta ~ normal(0,1);
  sigma ~ normal(0, 1);
  tau ~ normal(0, 1); 
  eps ~ normal(0, 1); 
  zeta ~ normal(0, 1);
  
  // calculate likelihood
  for (i in 1:C) {
    int lid = colony_land[i]; // landscape the colony belongs to
    int start_tr = traps_start[lid]; // start index for traps in that landscape
    int num_tr = traps_n[lid];  // number of traps in that landscape
    vector[num_tr] lambda_row; // place holder for visitation intensity of colony at each trap in the landscape
    int y_seg[num_tr]; // place holder for yobs values for colony i
    
    // PRIOR ON DELTA_X AND DELTA_Y
    // 5.841 is 99% crit value for student's t with df = 3...
    // we expect 99% of colonies to be within 2.5 km of site centroid
    // delta_x[i] ~ student_t(3, site_centroids[lid, 1], 2.5/5.841);
    // delta_y[i] ~ student_t(3, site_centroids[lid, 2], 2.5/5.841);

    for (t in 1:num_tr) {
      // fill y_seg with observations for that colony
      y_seg[t] = y_flat[ y_start[i] + t - 1 ];
      
      // compute lambda for each trap in that landscape
      int k = start_tr + t - 1;                 // global trap index
      real dis = sqrt( square(delta_x[i] - trap_pos[k,1]) +
                       square(delta_y[i] - trap_pos[k,2]) );
                       
      // penalize distances that are outside Rmax
      // a bit like a soft indicator function, to maintain differentiability for HMC
      target += -penalty * log1p_exp((dis-Rmax)*steepness);
      
      // calculate visitation intensity
      lambda_row[t] = -0.5*(dis / rho)^2 + beta*sample_effort[k] + zeta_scale[i] + eps_scale[k];
    }

    // compute multinomial probabilities and add to target likelihood
    vector[num_tr] multi_probs = softmax(lambda_row);
    y_seg ~ multinomial(multi_probs);
  }
}

