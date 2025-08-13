functions {
  vector uniform_point_rng(real lowerbound,
                           real upperbound_x,
                           real upperbound_y) {
    vector[2] x;
    
    x[1] = uniform_rng(lowerbound, upperbound_x);
    x[2] = uniform_rng(lowerbound, upperbound_y);
    
    return x;
  }
  
}

data {
}

transformed data {
  int<lower=0> C = 10000; // Number of colonies
  int<lower=0> K = 25;      // Number of traps
  int<lower=1> L = 6; // Number of landscapes
  int<lower=0> landscapesize = 1500; // Size of each landscape
  int<lower=0> trapgridsize = 300; // Size of each trapping grid
  real gridsize = sqrt(K); // Number of rows/columns in trap grid
  real stepsize = trapgridsize/(gridsize-1); // Distance between traps in grid
  
  real<lower=0> lowerbound = 0; // lower limit for colony locations
  real<lower=0> upperbound_x = 3*landscapesize;   // upper limit for colony locations (x direction)
  real<lower=0> upperbound_y = (L/3)*landscapesize;   // upper limit for colony locations (y direction)
  
  real alpha = log(1);
  real<lower=0> rho = 50;
  real theta = 0.5;
}

generated quantities {
  array[C] vector[2] delta;
  array[K*L] vector[2] trap_pos;
  vector[K*L] fq;
  
  array[C, K*L] int<lower=0> yobs;
  array[C, K*L] real<lower=0> sq_dist;
  
  for (i in 1:C)
    delta[i] = uniform_point_rng(lowerbound, upperbound_x, upperbound_y);
  
  for (l in 1:L) { // iterate over landscapes
    for (k in 1:K) { // iterate over traps
    
      //define landscape position
      int lancol = ((l-1) % 3);
      real lanrow = ceil(l/3.0) - 1;
      
      // define trap position
      int trapnum = k + (l-1)*K;
      int trapcol = (k-1) % to_int(gridsize);
      real traprow = ceil(k/gridsize) - 1;
      
      // get trap coordinates
      trap_pos[trapnum][1] = lancol*landscapesize + (landscapesize - trapgridsize)/2 + (trapcol)*stepsize;
      trap_pos[trapnum][2] = lanrow*landscapesize + (landscapesize - trapgridsize)/2 + (traprow)*stepsize;
      fq[trapnum] = normal_rng(0,1);

      
      for (i in 1:C) { // iterate over colonies
        real d = squared_distance(trap_pos[trapnum], delta[i]);
        sq_dist[i, trapnum] = d;
        yobs[i, trapnum] = poisson_log_rng(alpha - d / square(rho) + theta*fq[trapnum]);
      }
    }
    
  }
  
}
