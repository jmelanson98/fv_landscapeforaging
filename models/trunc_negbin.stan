// Fit truncated poisson!
// from https://freerangestats.info/blog/2018/03/20/truncated-poisson
//

data {
  int num_sites; // number of sites surveyed
  int num_colonies; // the total number of colonies at all sites
  int lower_limit; // the lower limit of detections, e.g., 1 individual per colony
  int <lower = lower_limit> sib_counts[num_colonies]; // the number of siblings observed for each colony
  int <lower = 1, upper = num_sites> site_ids[num_colonies]; // the site id for each colony
  real alpha_mu;
  real alpha_sigma;
}

parameters {
  vector<lower = 0>[num_sites] alpha;
  vector<lower = 0>[num_sites] beta;
}

model {
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(0,1);
  
  for(i in 1:num_colonies){
    sib_counts[i] ~ neg_binomial(alpha[site_ids[i]], beta[site_ids[i]]) T[lower_limit, ];
  }
}

