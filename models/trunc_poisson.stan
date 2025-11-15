// Fit truncated poisson!
// from https://freerangestats.info/blog/2018/03/20/truncated-poisson
//

data {
  int num_sites; // number of sites surveyed
  int num_colonies; // the total number of colonies at all sites
  int lower_limit; // the lower limit of detections, e.g., 1 individual per colony
  int <lower = lower_limit> sib_counts[num_colonies]; // the number of siblings observed for each colony
  int <lower = 1, upper = num_sites> site_ids[num_colonies]; // the site id for each colony
  real lambda_mu;
  real lambda_sigma;
}

parameters {
  vector<lower = 0>[num_sites] lambda; 
}

model {
  lambda ~ normal(lambda_mu, lambda_sigma);
  
  for(i in 1:num_colonies){
    sib_counts[i] ~ poisson(lambda[site_ids[i]]) T[lower_limit, ];
  }
}



generated quantities {
  // probability of observing a zero
  real pzero = exp(poisson_lpmf(0 | lambda));
  
  // total_colonies = num_colonies / P(observed). 
  // P(observed) = 1 - p_zero
  real total_colonies = num_colonies / (1 - pzero);
}
