// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=0> C;                // total number of colonies
  int<lower=0> K;                // total number of traps
  int<lower=0> O;                // total number of observations (observation for a colony at a trap, in a timepoint)
  int<lower=0> CT;               // number of colonies * number of timepoints
  int CTstarts[CT];              // start index for each colony-trap combo
  int CTlengths[CT];             // number of obs per colony-trap combo
  matrix[O,2] trap_pos;              // trap coordinates
  int colony_id[O];              // colony id for each observation
  int trap_id[O];                 // trap id for each observation
  int landscape_id[O];            // landscape id for each observation
  vector[O] fq;                      // floral quality for each observation
  int yobs[O];                   // number of individuals observed
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
  real<lower=0> beta;
  real<lower=0> sigma;
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
  beta ~ normal(0,1);
  sigma ~ normal(0, 1);
  eps ~ normal(0, 1);
  
  // calculate likelihood
  for (n in 1:CT){
    int start = CTstarts[n]; // index start
    int length = CTlengths[n]; // index length
    
    // get data subsets
    int y_it[length] = yobs[start : start+length-1]; // yobs for colony i, timepoint t set
    vector[length] fq_it = fq[start:start+length-1]; // yobs for colony i, timepoint t set
    int k[length] = trap_id[start:start+length-1]; // trap indices
    int i = colony_id[start]; // colony id
    matrix[length,2] trap_t = trap_pos[start:start+length-1,];
    
    // compute lambda for each trap in that landscape
    vector[length] dis = sqrt( square(delta_x[i] - trap_t[,1]) +
                       square(delta_y[i] - trap_t[,2]) );
    vector[length] lambda_it = -0.5*(dis / rho)^2 + beta*fq_it + eps_scale[k];
    
    // compute multinomial probabilities and add to target likelihood
    vector[length] multi_probs = softmax(lambda_it);
    y_it ~ multinomial(multi_probs);
                       
    // penalize distances that are outside Rmax
    // a bit like a soft indicator function, to maintain differentiability for HMC
    target += to_vector(y_it) .* (-penalty * log1p_exp((dis-Rmax)*steepness));
  }
}
