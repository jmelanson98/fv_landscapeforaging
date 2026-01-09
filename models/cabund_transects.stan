data {
  int<lower=1> O;         // # observations
  int<lower=1> T;         // # transects
  vector[O] idw;            // trap idw
  vector[O] iji;            // trap iji
  vector[O] fq;             // floral quality
  vector[O] year;           // 0/1 = 2022/2023
  vector[O] effort;          // sampling effort
  int transect_id[O];
  int yobs[O];            // number of colonies observed
}

parameters {
  real alpha; //global intercept
  real beta_idw;
  real beta_iji;
  real beta_fq;
  real beta_yr;
  real<lower=0> phi;
  real<lower=0> sigma;
  vector[T] eps;
}

transformed parameters {
  // get eps_scale
  vector[T] eps_scale = eps*sqrt(sigma);
  
  vector[O] mu;
  for (n in 1:O) {
    mu[n] = exp(
      alpha
      + beta_idw * idw[n]
      + beta_iji * iji[n]
      + beta_fq * fq[n]
      + beta_yr * year[n]
      + eps_scale[transect_id[n]]
      + log(effort[n])
    );
  }
}


model {
  // set priors
  alpha ~ normal(2.5,1);
  beta_idw ~ normal(0,1);
  beta_iji ~ normal(0,1);
  beta_fq ~ normal(0,3);
  beta_yr ~ normal(0,3);
  phi ~ exponential(1);
  sigma ~ normal(0, 1);
  eps ~ normal(0, 1);
  
  
  // likelihood
  target += neg_binomial_2_lpmf(yobs | mu, phi);
}



