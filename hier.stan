// The user-defined functions
functions {
  // Function that encapsulates the dynamics of the system
  vector dx_dt(real t, vector x, real m, real l, real mu, real g) { 
    vector[2] dxdt;
    
    dxdt[1] = x[2];
    dxdt[2] = (-mu/(m*l))*x[2] -  (g/l)*x[1];
    return dxdt;
  }
  
  // Function that solves the ODE and returns the position array
  real[] simulate_x(int N, real[] t, vector s0, real m,
                    real l, real mu, real g){
    vector[2] z[N-1] = ode_rk45(dx_dt, s0, 0.0, t[2:N], m, l, mu, g);
    array[1] real x0;
    x0[1] = s0[1];
    real x_sim[N] = append_array(x0, z[,1]);
    
    return x_sim;
  }
}

// The problem data
data {
  int<lower=1> N;
  matrix[N,2] x;
  array[N,2] real ts;
  matrix[9, 2] pr_params; // Matrix for prior parameters
}

// The problem parameters
// <lower=0.000001> is to prevent division by zero issues
parameters {
  vector[2] state01;
  vector[2] state02;
  real<lower=0> sigma;
  real<lower=0.000001> m0; // mean mass
  vector<lower=0.000001>[2] m; // masses for both nuts 
  real<lower=0.0001> l; 
  real<lower=0> mu0; // mean air-resistance
  vector<lower=0>[2] mu; // air-resistance for both nuts
  real<lower=0> g;
}

// The ODE-transformation from params to simulated data
transformed parameters {
  real x_sim1[N] = simulate_x(N, ts[,1], state01, m[1], l, mu[1], g);
  real x_sim2[N] = simulate_x(N, ts[,2], state02, m[2], l, mu[2], g);
}

// The actual model
model {
  // The priors
  m0 ~ normal(pr_params[1,1], pr_params[9,1]);
  l ~ normal(pr_params[2,1], pr_params[2,2]);
  mu0 ~ uniform(pr_params[3,1], pr_params[3,2]);
  g ~ normal(pr_params[4,1], pr_params[4,2]);
  state01[1] ~ normal(pr_params[5,1], pr_params[5,2]); // Position of nut 1
  state01[2] ~ normal(pr_params[7,1], pr_params[7,2]); // Velocity of nut 1
  state02[1] ~ normal(pr_params[6,1], pr_params[6,2]); // Position of nut 2
  state02[2] ~ normal(pr_params[7,1], pr_params[7,2]); // Velocity of nut 2
  sigma ~ normal(pr_params[8,1], pr_params[8,2]);
  
  // Hierarchical part
  for (n in 1:2) {
    m[n] ~ normal(m0, pr_params[1,2]);
    mu[n] ~ normal(mu0, 2*mu0);
  }
  
  // The likelihood model
  for (t in 1:N) {
    x[t,1] ~ normal(x_sim1[t], sigma);
    x[t,2] ~ normal(x_sim2[t], sigma);
  }
  
}

// Computing the needed quantities for LOO-CV and other statistics
generated quantities {
  // Posterior draws
  real x_draws1[N] = simulate_x(N, ts[,1], state01, m[1], l, mu[1], g);
  real x_draws2[N] = simulate_x(N, ts[,2], state02, m[2], l, mu[2], g);
  
  // Log-likelihoods
  matrix[N,2] log_lik;
  for (t in 1:N){
    log_lik[t,1] = normal_lpdf(x[t,1] | x_sim1[t], sigma);
    log_lik[t,2] = normal_lpdf(x[t,2] | x_sim2[t], sigma);
  }
}

