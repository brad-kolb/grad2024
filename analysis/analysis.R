#### packages ####
library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(cowplot)

#### compile stan program that specifies the model ####
model <- cmdstan_model(here("models", "model.stan"))

#### compile stan program that specifies the generated quantities ####
model_gq <- cmdstan_model(here("models", "generated_quantities.stan"))

#### import data ####
data <- read_csv(file = here("data", "clean_data.csv")) %>% 
  as_tibble() %>% 
  select(J, K, n_c, r_c, n_t, r_t) %>% 
  filter(K == 4) %>% 
  mutate(j = row_number()) %>% 
  pivot_longer(cols = starts_with(c("n_", "r_")), 
               names_to = c(".value", "arm"), 
               names_sep = "_") %>% 
  mutate(x = ifelse(arm == "c", 0, 1),
         l = row_number()) %>%
  select(K, j, l, x, n, r)

dat <- with(data,
            list(J = length(unique(j)),
                 L = length(l),
                 jj = j,
                 ll = l,
                 x = x,
                 n = n,
                 y = r,
                 estimate_posterior = 1, 
                 priors = 1))

fit <- model$sample(data = dat, seed = 1, chains = 4, 
                    parallel_chains = 4, save_warmup = TRUE, refresh = 1000)

#### evaluate convergence diagnostics ####

# rhat should not in general exceed 1.01 for each parameter
map(c("rho", "sigma", "mu", "tau"), 
    function(params) {
      posterior::rhat(extract_variable_matrix(fit, variable = params))
    }
)

# bulk ESS should in general exceed 100 for each parameter
map(c("rho", "sigma", "mu", "tau"), 
    function(params) {
      posterior::ess_bulk(extract_variable_matrix(fit, variable = params))
    }
)

#### obtain generated quantities ####
fit_gq <- model_gq$generate_quantities(fit, data = dat, seed = 123)

#### get estimates ####

# what are my chances without treatment?
fit_gq$summary(
  variables = "control_marg",
  median, mad,
  extra_quantiles = ~posterior::quantile2(., probs = c(0.025, .975))
)

# what are my chances with treatment?
fit_gq$summary(
  variables = "treatment_marg",
  median, mad,
  extra_quantiles = ~posterior::quantile2(., probs = c(0.025, .975))
)

# how certain are you that my chances improve with treatment
fit_gq$summary(
  variables = c("rrr_marg", "arr_marg"),
  median, mad,
  extra_quantiles = ~posterior::quantile2(., probs = c(0.025, .975)),
  pr_pos = ~ mean(. > 0)
  )

#### make plots ####

# average risk
control <- fit_gq %>%
  spread_draws(control_marg) %>%
  rename(value = control_marg) %>% 
  mutate(type = "Medical")

treatment <- fit_gq %>% 
  spread_draws(treatment_marg) %>% 
  rename(value = treatment_marg) %>% 
  mutate(type = "Surgical")

draws <- rbind(control, treatment)

risk_intervals <- draws %>%
  ggplot(aes(x = type, y = value)) +
  stat_pointinterval(.width = .95) +
  labs(title = "Average risk of poor outcome, all trials",
       y = "Risk, %",
       x = NULL) +
  scale_y_continuous(limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  cowplot::theme_half_open()

ggsave(here("analysis", "risk_intervals.png"), 
       risk_intervals, 
       width = 8, height = 6, dpi = 300, bg = "white")

risk_densities <- draws %>%
  ggplot(aes(x = type, y = value)) +
  stat_slabinterval(.width = .95, side = "left") +
  labs(title = "Average risk of poor outcome, all trials",
       y = "Risk, %",
       x = NULL) +
  scale_y_continuous(limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  cowplot::theme_half_open()

ggsave(here("analysis", "risk_densities.png"), 
       risk_densities, 
       width = 8, height = 6, dpi = 300, bg = "white")

# risk reduction
arr <- fit_gq %>% 
  spread_draws(arr_marg) %>% 
  rename(value = arr_marg) %>% 
  mutate(type = "Absolute")

rrr <- fit_gq %>% 
  spread_draws(rrr_marg) %>% 
  rename(value = rrr_marg) %>% 
  mutate(type = "Relative")

draws <- rbind(arr, rrr)

reductions_intervals <- draws %>%
  ggplot(aes(x = type, y = value)) +
  stat_slabinterval(.width = .95, side = "left") +
  labs(title = "Average risk reduction, absolute vs relative",
       y = "Reduction, %",
       x = NULL) +
  cowplot::theme_half_open() +
  scale_y_continuous(limits = c(-.1,.4), labels = c(-10, 0, 10, 20, 30, 40)) +
  guides(fill = "none")

ggsave(here("analysis", "reductions_intervals.png"), 
       reductions_intervals, 
       width = 8, height = 6, dpi = 300, bg = "white")

reductions_probabilities <- draws %>%
  ggplot(aes(x = type, y = value)) +
  stat_slab(aes(fill = after_stat(y > 0)), side = "left") +
  labs(title = "Average risk reduction, absolute vs relative",
       y = "Reduction, %",
       x = NULL) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  cowplot::theme_half_open() +
  scale_y_continuous(limits = c(-.1,.4), labels = c(-10, 0, 10, 20, 30, 40)) +
  guides(fill = "none")

ggsave(here("analysis", "reductions_probabilities.png"), 
       reductions_probabilities, 
       width = 8, height = 6, dpi = 300, bg = "white")
  