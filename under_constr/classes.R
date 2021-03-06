# library(Rcpp)
# devtools::clean_dll()
# Rcpp::compileAttributes()
# devtools::load_all()
library(sappp)
library(tidyverse)

# if (is.list(x) & is.null(x$a)) print("cool")


k <- runif(8)
r <- runif(8)
so <- make_sim_obj(
    .n_pops = 8, .n_patches = 1,
    .harvest_periods = 1e6, .harvest_offsets = 0,
    wasp_density_0 = 0,
    # aphid_surv_juv = as.list(sappp::populations$surv_juv$high),
    # aphid_surv_adult = lapply(r,
    #                           function(x) sappp::populations$surv_adult$high * x),
    K = as.list(sappp::populations$K * k),
    disp_aphid = list(0.75), harvest_surv = 0.1)

sim <- so$simulate(200)

sim_df <- sim$flatten() %>% 
    as_data_frame %>%
    rename(patch = V1, line = V2, t = V3, aphids = V4, parasit = V5, 
           mummies = V6, wasps = V7) %>% 
    mutate(patch = factor(as.integer(patch)), 
           line = factor(as.integer(line)),
           t = as.integer(t))
# source(".Rprofile")
sim_df %>% 
    # filter(t <= 50) %>%
    # filter(patch == 1) %>% 
    ggplot(aes(t, aphids, color = line)) +
    theme_classic() + 
    theme(strip.background = element_blank()) +
    facet_wrap(~ patch, nrow = 2) +
    geom_line() +
    geom_line(data = group_by(sim_df, patch, t) %>% summarize(wasps = sum(wasps)),
              aes(y = max(sim_df$aphids) * wasps / max(wasps)), 
              color = 'black', linetype = 2)



set.seed(2)
so <- make_sim_obj(
    .n_pops = 10, .n_patches = 8,
    .harvest_periods = 7,
    .harvest_offsets = c(0, 3),
    aphid_repro = list(.sd_pops = 0, .sd_patches = 0.5),
    aphid_surv_juv = list(.sd_pops = 0.1),
    aphid_surv_adult = list(.sd_pops = 0.1),
    aphid_repro = list(.sd_pops = 0.01),
    K = list(.sd_pops = 0.01),
    # K = NULL,
    disp_aphid = list(.sd_pops = 0.25, .mean = 0.5),
    attack_surv = list(.sd_pops = 0.5),
    no_error = TRUE,
    harvest_surv = 0.5,
    wasp_density_0 = 1)
sim <- so$simulate(120)

sim_df <- sim$flatten() %>% 
    as_data_frame %>%
    rename(patch = V1, line = V2, t = V3, aphids = V4, parasit = V5, 
           mummies = V6, wasps = V7) %>% 
    mutate(patch = factor(as.integer(patch)), 
           line = factor(as.integer(line)),
           t = as.integer(t))

source(".Rprofile")
sim_df %>% 
    # filter(t <= 50) %>%
    # filter(patch == 1) %>% 
    ggplot(aes(t, aphids, color = line)) +
    theme_classic() + 
    theme(strip.background = element_blank()) +
    facet_wrap(~ patch, nrow = 2) +
    geom_line() +
    geom_line(data = group_by(sim_df, patch, t) %>% summarize(wasps = sum(wasps)),
              aes(y = max(sim_df$aphids) * wasps / max(wasps)), 
              color = 'black', linetype = 2)


