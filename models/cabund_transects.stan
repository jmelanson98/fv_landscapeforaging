data {
  int<lower=1> T;         // # transects
  vector[T] idw;            // trap idw
  vector[T] iji;            // trap iji
  vector[T] fq;             // floral quality
  vector[T] effort;          // sampling effort
  int yobs[T];            // number of colonies observed
}

parameters {
  real alpha; //global intercept
  real beta_idw;
  real beta_iji;
  real beta_fq;
  real<lower=0> beta_eff;
  real<lower=0> phi;
}

transformed parameters {
  vector[T] mu;
  mu = exp(alpha + beta_idw*idw + beta_iji*iji + beta_fq*fq + beta_eff*effort);  // log link
}


model {
  // set priors
  alpha ~ normal(2.5,1);
  beta_idw ~ normal(0,1);
  beta_iji ~ normal(0,1);
  beta_fq ~ normal(0,3);
  beta_eff ~ normal(0,1);
  phi ~ exponential(1);
  
  // likelihood
  target += neg_binomial_2_lpmf(yobs | mu, phi);
}



