//Modified from Pope & Jha, 2017, Cons. Genetics

data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
matrix[K, 2] trap_pos; // trap coordinates 
int y[C, K]; // counts of bees in traps
real lowerbound; // uniform prior on colony location 
real upper_y; // uniform prior on colony location
real upper_x; // uniform prior on colony location
vector[K] floral; // floral quality at traps
}

parameters { // see text for definitions
real<lower=0> rho; 
real<lower=0> sigma;
real<lower=0> tau; 
real theta;
vector[K] eps; 
vector[C] zeta;
array[C] real<lower=lowerbound, upper=upper_x> delta_x;
array[C] real<lower=lowerbound, upper=upper_y> delta_y;


}


transformed parameters {
  real<lower=0> tau_sqrt = sqrt(tau);
  real<lower=0> sigma_sqrt = sqrt(sigma); 
  vector[C] zeta_scale = zeta * tau_sqrt;
  vector[K] eps_scale = eps*sigma_sqrt;
}
  
  
model {
  
// temporary declarations
matrix[C,K] dis; //distance from colony C to trap K
matrix[C,K] lambda; //rate of captures for colony C at trap K
vector[K] multi_probs; //multinomial probability of capture from colony C at trap K

// priors
sigma ~ normal(0, 1);
tau ~ normal(0, 1); 
rho ~ lognormal(3.5, 0.5);
theta ~ normal(0, 1);

// random effects for traps
eps ~ normal(0, 1); 
zeta ~ normal(0, 1);

// calculate intensity
    for(i in 1:C){
      for(k in 1:K){
        dis[i, k] = sqrt(square(delta_x[i] - trap_pos[k,1]) + square(delta_y[i] - trap_pos[k,2]));
        lambda[i, k] = (-dis[i,k]/rho) + theta*floral[k] + zeta_scale[i] + eps_scale[k];
      }
      // compute  multinomial probabilities for colony i
      multi_probs = softmax(lambda[i,]');
    
      // add to target likelihood (trap counts for colony i)
      y[i,] ~ multinomial(multi_probs);
    }
}

// generated quantities {
//   vector[C] colony_dist;        // Declare estimated colony foraging distance
//   colony_dist = rep_vector(0, C);  // Initialize colony distance to 0
//   
//   // start anonymous scope
//   {
//     matrix[C,K] dis; //distance from colony C to trap K
//     matrix[C,K] lambda; //rate of captures for colony C at trap K
//     vector[C] V;     // Declare local variable for each colony's total visitation
//     real alpha = 1e-12;  // set smaaaaall number to prevent division by 0
//     
//     // Recompute lambda and dis
//     for(k in 1:K){
//       for(i in 1:C){
//         dis[i, k] = sqrt(square(delta[i, 1] - trap[k,1]) + square(delta[i, 2] - trap[k,2]));
//         lambda[i, k] = dis[i,k]/(-rho*exp(theta*floral[k])) + mu + zeta_scale[i] + eps_scale[k];
//       } 
//     }
//     
//     // Compute V for normalization
//     for (i in 1:C){
//       V[i] = sum(exp(lambda[i,]));
//     }
//     
//     
//     // compute colony_dist to be saved outside the anonymous scope
//     for (k in 1:K){
//       colony_dist = colony_dist + (dis[,k] .* exp(lambda[,k]) ./ (V + alpha));
//     }
//   }
// }
// 
