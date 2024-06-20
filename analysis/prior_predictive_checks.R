#### packages ####
library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(cowplot)

#### options ####
theme_set(new = theme_minimal() +
            theme(
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank()
            ))

##### compile stan program that specifies the model ####
model <- cmdstan_model(here("models", "model.stan"))

#### compile stan program that specifies the generated quantities ####
model_gq <- cmdstan_model(here("models", "generated_quantities.stan"))

# obtain predictors
predictors <- read_csv(file = here("data", "clean_data.csv")) %>% 
  as_tibble() %>% 
  select(J, K, n_c, r_c, n_t, r_t) %>% 
  pivot_longer(cols = starts_with(c("n_", "r_")), 
               names_to = c(".value", "arm"), 
               names_sep = "_") %>% 
  mutate(x = ifelse(arm == "c", 0, 1),
         l = row_number(),
         j = J) %>%
  select(K, j, l, x, n) %>% 
  filter(K == 1)
  
dat <- with(predictors,
            list(J = length(unique(j)),
                 L = length(l),
                 jj = j,
                 ll = l,
                 x = x,
                 n = n,
                 y = rep(0, length(l)),
                 estimate_posterior = 0, 
                 priors = 1))

pps1 <- model$sample(data = dat, seed = 123, chains = 4, 
                     parallel_chains = 4, save_warmup = TRUE, refresh = 1000)

pps1_gq <- model_gq$generate_quantities(pps1, data = dat, seed = 123)

dat$priors <- 0

pps0 <- model$sample(data = dat, seed = 123, chains = 4, 
                     parallel_chains = 4, save_warmup = TRUE, refresh = 1000)

pps0_gq <- model_gq$generate_quantities(pps0, data = dat, seed = 123)

# expected trial specific event probabilities for control and treatment arms
# implied by the prior distributions for the model parameters

bayesplot::mcmc_hist(pps0_gq$draws("E_y_tilde"))
bayesplot::mcmc_hist(pps0_gq$draws("E_rr_tilde"))
# bayesplot::mcmc_hist(pps0_gq$draws("E_y_tilde"))

pps1_gq$summary() %>% print(n = Inf)
pps0_gq$summary() %>% print(n = Inf)
