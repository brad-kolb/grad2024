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
generated quantities{
  // expected trial-specific event probability for treatment and control arms
  vector<lower=0,upper=1>[L] E_y_tilde;
  for (i in 1:L) {
    E_y_tilde[i] = 1 - inv_logit(phi[jj[i]] + theta[jj[i]] * x[i]);
  }
  // expected trial-specific relative risk 
  vector[J] E_arr_tilde;
  for (i in 1:J) {
    E_arr_tilde[i] = E_y_tilde[2 * i - 1] - E_y_tilde[2 * i];
  }
  vector[J] E_rrr_tilde;
  for (i in 1:J) {
    E_rrr_tilde[i] = (E_y_tilde[2 * i - 1] - E_y_tilde[2 * i]) / E_y_tilde[2 * i - 1];
  }
  // marginal effect of intervention across trials
  real control_marg = 1 - inv_logit(rho);
  real treatment_marg = 1 - inv_logit(rho + mu);
  real arr_marg = control_marg - treatment_marg;
  real rrr_marg = (control_marg - treatment_marg) / control_marg;
  // expected effect of intervention in the next hypothetical trial
  real E_control_next = 1 - inv_logit(normal_rng(rho, sigma));
  real E_treatment_next = 1 - inv_logit(normal_rng(rho, sigma) + normal_rng(mu, tau));
  real E_arr_next = E_control_next - E_treatment_next;
  real E_rrr_next = (E_control_next - E_treatment_next) / E_control_next;
}
