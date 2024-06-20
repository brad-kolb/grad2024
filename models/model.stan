data {
  int<lower=1> J; // total number of studies
  int<lower=1> L;  // total number of arms (each study has two arms)
  array[L] int<lower=1,upper=J> jj;  // study ID
  array[L] int<lower=1,upper=L> ll; // arm ID
  array[L] int<lower=0,upper=1> x;  // intervention covariate 
  array[L] int<lower=0> y;  // number of events in each arm
  array[L] int<lower=0> n;  // number of observations in each arm
  int<lower=0> estimate_posterior;  // switch for estimating posterior vs running prior predictive simulation
  int<lower=0> priors;  // switch for checking sensitivity of posterior to alternative specification for priors
}
parameters {
  real rho;  // population mean baseline log odds
  real<lower=0> sigma;  // population sd of baseline log odds
  vector<offset=rho, multiplier=sigma>[J] phi;  // per study baseline log odds
  real mu;  // population mean treatment effect (difference in log odds of outcome at baseline and with treatment)
  real<lower=0> tau;  // population sd of treatment effects
  vector<offset=mu, multiplier=tau>[J] theta;  // per study treatment effect (log odds difference)
}
model {
  if (estimate_posterior == 1) {
    // linear model
    vector[L] q;
    for (i in 1:L) {
      q[i] = phi[jj[i]] + theta[jj[i]] * x[i];
    }
    // likelihood
    y[1:L] ~ binomial(n[1:L], inv_logit(q[1:L]));
  }
  // hyperpriors
  phi[1:J] ~ normal(rho, sigma);
  theta[1:J] ~ normal(mu, tau);
  // priors
  if (priors == 1) {
  rho ~ normal(-1, 1);
  sigma ~ normal(0, 1);
  mu ~ normal(0, 1);
  tau ~ normal(0, 1);
  } else { 
    mu ~ std_normal(); 
    rho ~ std_normal();
    sigma ~ std_normal();
    tau ~ std_normal();
  }
}
