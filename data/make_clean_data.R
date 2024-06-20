# load packages -----------
library(here)
library(tidyverse)

# import raw data --------
# this data was manually entered into a spreadsheet by me
data <- read_csv(file = here("data", "raw_data.csv"))

# define functions to process raw data ------------
make_clean_data <- function(data, metric) {
  df <- tibble(
    J = data %>% # trials
      filter(treatment_id == 1) %>% 
      .$trial_id, 
    K = data %>% # trial types
      filter(treatment_id == 1) %>% 
      .$group_id, 
    n_c = data %>% # number of cases, control
      filter(treatment_id == 0) %>% 
      .$tot_actual,
    r_c = data %>% # number of successes, control
      filter(treatment_id == 0) %>%
      .[[metric]], # dynamically accessing the column
    n_t = data %>% # number of cases, treatment
      filter(treatment_id == 1) %>% 
      .$tot_actual,
    r_t = data %>% # number of successes, treatment
      filter(treatment_id == 1) %>%
      .[[metric]] # dynamically accessing the column
  ) %>% 
    mutate(
      y = log(r_t / (n_t - r_t)) - log(r_c / (n_c - r_c)), # log odds ratio
      se = ifelse(r_t == 0 | r_c == 0, # approximate standard error
                  NA_real_,
                  sqrt(1/r_t + 1/(n_t - r_t) + 1/r_c + 1/(n_c - r_c))
      )
    )
  return(df)
}

# process raw data and save output as csv file ---------
df <- make_clean_data(data, "independent")

df %>% write_csv(
  here("data", "clean_data.csv")
)