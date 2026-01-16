// Lineage turnover mixture model!
// J Melanson
// started January 14, 2026

data {
  int<lower=0> C; // number of lineages (2022 colonies)
  int y[C]; // number of queens observed for each lineage in 2023
}

parameters {
  real<lower=0, upper=1> theta;
  real<lower=0> lambda;
}


model {
  theta ~ beta(9,1);
  lambda ~ normal(0, 1);
  
  for (c in 1:C){
    if (y[c] == 0) {
      target += log_sum_exp(log(theta),
                            log1m(theta)
                              + poisson_lpmf(y[c] | lambda));
    } else {
      target += log1m(theta)
                  + poisson_lpmf(y[c] | lambda  );
    }
  }
}

