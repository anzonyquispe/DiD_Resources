# Load libraries and set baseline parameters
library(tidyverse)
library(lfe)
library(fastDummies)
library(ggthemes)
library(did)
theme_set(theme_clean() + theme(plot.background = element_blank()))
#----------------------------------------------------------------------------
iseed  = 20201221
nrep <- 100  
true_mu <- 1
set.seed(iseed)

#----------------------------------------------------------------------------
## Generate data - treated cohorts consist of 250 obs each, with the treatment effect still = true_mu on average
make_data <- function(nobs = 1000, 
                      nstates = 40) {
  
  # unit fixed effects (unobservd heterogeneity)
  unit <- tibble(
    unit = 1:nobs,
    # generate state
    state = sample(1:nstates, nobs, replace = TRUE),
    unit_fe = rnorm(nobs, state/5, 1),
    # generate instantaneous treatment effect
    #mu = rnorm(nobs, true_mu, 0.2)
    mu = true_mu
  )
  
  # year fixed effects (first part)
  year <- tibble(
    year = 1980:2010,
    year_fe = rnorm(length(year), 0, 1)
  )
  
  # Put the states into treatment groups
  treat_taus <- tibble(
    # sample the states randomly
    state = sample(1:nstates, nstates, replace = FALSE),
    # place the randomly sampled states into four treatment groups G_g
    cohort_year = sort(rep(c(1986, 1992, 1998, 2004), 10))
  )
  
  # make main dataset
  # full interaction of unit X year 
  expand_grid(unit = 1:nobs, year = 1980:2010) %>% 
    left_join(., unit) %>% 
    left_join(., year) %>% 
    left_join(., treat_taus) %>% 
    # make error term and get treatment indicators and treatment effects
    # Also get cohort specific trends (modify time FE)
    mutate(error = rnorm(nobs*31, 0, 1),
           treat = ifelse(year >= cohort_year, 1, 0),
           tau = ifelse(treat == 1, mu, 0),
           year_fe = year_fe + 0.1*(year - cohort_year)
    ) %>% 
    # calculate cumulative treatment effects
    group_by(unit) %>% 
    mutate(tau_cum = cumsum(tau)) %>% 
    ungroup() %>% 
    # calculate the dep variable
    mutate(dep_var = (2010 - cohort_year) + unit_fe + year_fe + tau_cum + error)
  
}
#----------------------------------------------------------------------------
# make data
data <- make_data()

# plot
plot1 <- data %>% 
  ggplot(aes(x = year, y = dep_var, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = data %>% 
              group_by(cohort_year, year) %>% 
              summarize(dep_var = mean(dep_var)),
            aes(x = year, y = dep_var, group = factor(cohort_year),
                color = factor(cohort_year)),
            size = 2) + 
  labs(x = "", y = "Value", color = "Treatment group   ") + 
  geom_vline(xintercept = 1986, color = '#E41A1C', size = 2) + 
  geom_vline(xintercept = 1992, color = '#377EB8', size = 2) + 
  geom_vline(xintercept = 1998, color = '#4DAF4A', size = 2) + 
  geom_vline(xintercept = 2004, color = '#984EA3', size = 2) + 
  scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom',
        #legend.title = element_blank(), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))  +
  ggtitle("One draw of the DGP with homogeneous effects across cohorts \n and with all groups being eventually treated")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot1

#ggsave("plot_dgp1.png", plot1, width = 10, height = 5, dpi = 100)





# function to run ES DID
# variables we will use
keepvars <- c("`rel_year_-5`",  "`rel_year_-4`",  "`rel_year_-3`",  "`rel_year_-2`",
              "rel_year_0", "rel_year_1", "rel_year_2", "rel_year_3", "rel_year_4", "rel_year_5")

run_ES_DiD <- function(...) {
  
  # resimulate the data
  data <- make_data()
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year) %>% 
    dummy_cols(select_columns = "rel_year") %>% 
    # generate pre and post dummies
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           Post = ifelse(rel_year > 5, 1, 0))
  
  # estimate the model
  mod <- lfe::felm(dep_var ~ Pre + `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + `rel_year_-2` + 
                     `rel_year_0` + `rel_year_1` + `rel_year_2` + `rel_year_3` + `rel_year_4` + 
                     `rel_year_5` + Post | unit + year | 0 | state, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  mod2 <- tibble(
    estimate = mod$coefficients,
    term1 = rownames(mod$coefficients)
  )
  
  es <-
    mod2 %>% 
    filter(term1 %in% keepvars) %>% 
    mutate(t = c(-5:-2, 0:5)) %>% 
    select(t, estimate)
  es
}
run_ES_DiD2 <- function(...) {
  
  # resimulate the data
  data <- make_data()
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year) %>% 
    dummy_cols(select_columns = "rel_year") %>% 
    # generate pre and post dummies
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           Post = ifelse(rel_year > 5, 1, 0))
  
  # estimate the model
  mod <- feols(dep_var ~ Pre + `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + `rel_year_-2` + 
                     `rel_year_0` + `rel_year_1` + `rel_year_2` + `rel_year_3` + `rel_year_4` + 
                     `rel_year_5` + Post | unit + state^year, data = data, cluster = ~state)
  
  # grab the obs we need
  mod2 <- tibble(
    estimate = mod$coefficients,
    term1 = rownames(mod$coefficients)
  )
  
  es <-
    mod2 %>% 
    filter(term1 %in% keepvars) %>% 
    mutate(t = c(-5:-2, 0:5)) %>% 
    select(t, estimate)
  es
}
data_classical <- map_dfr(1:nrep, run_ES_DiD)
# data_classical2 <- map_dfr(1:3, run_ES_DiD2)

colors <- c("True Effect" = "red", "Estimated Effect" = "blue")

ES_plot_classical <- data_classical %>% 
  group_by(t) %>% 
  summarize(avg = mean(estimate),
            sd = sd(estimate),
            lower.ci = avg - 1.96*sd,
            upper.ci = avg + 1.96*sd) %>% 
  bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
  ggplot(aes(x = t, y = avg)) + 
  #geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
  geom_point(color = 'blue', size = 3) + 
  geom_line(aes(color = 'Estimated Effect'), size = 1) + 
  geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle("TWFE event-study regression with binned end-points")+
  scale_color_manual(values = colors) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom", 
        legend.title = element_blank())

ES_plot_classical



# function to run ES DID
run_ES_DiD_sat <- function(...) {
  
  # resimulate the data
  data <- make_data()
  
  # make dummy columns
  data <- data %>% 
    # make relative year indicator
    mutate(rel_year = year - cohort_year)
  
  # get the minimum relative year - we need this to reindex
  min_year <- min(data$rel_year)
  
  # reindex the relative years
  data <- data %>% 
    mutate(rel_year = rel_year - min_year) %>% 
    dummy_cols(select_columns = "rel_year")
  
  # make regression formula 
  indics <- paste("rel_year", (1:max(data$rel_year))[-(-1 - min_year)], sep = "_", collapse = " + ")
  keepvars <- paste("rel_year", c(-5:-2, 0:5) - min_year, sep = "_")  
  formula <- as.formula(paste("dep_var ~", indics, "| unit + year | 0 | state"))
  
  # run mod
  mod <- felm(formula, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  mod2 <- tibble(
    estimate = mod$coefficients,
    term1 = rownames(mod$coefficients)
  )
  
  es <-
    mod2 %>% 
    filter(term1 %in% keepvars) %>% 
    mutate(t = c(-5:-2, 0:5)) %>% 
    select(t, estimate)
  es
}

data_sat <- map_dfr(1:nrep, run_ES_DiD_sat)

ES_plot_sat <- data_sat %>% 
  group_by(t) %>% 
  summarize(avg = mean(estimate),
            sd = sd(estimate),
            lower.ci = avg - 1.96*sd,
            upper.ci = avg + 1.96*sd) %>% 
  bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
  ggplot(aes(x = t, y = avg)) + 
  #geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
  geom_point(color = 'blue', size = 3) + 
  geom_line(aes(color = 'Estimated Effect'), size = 1) + 
  geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle("TWFE event-study regression with 'all' leads and lags")+
  scale_color_manual(values = colors) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom", 
        legend.title = element_blank())


ES_plot_sat

#ggsave("es_plot_sat.png", ES_plot_sat, width = 10, height = 5, dpi = 150)
#ggsave("es_plot_classical.png", ES_plot_classical, width = 10, height = 5, dpi = 150)\


# function to run ES DID
run_CS <- function(...) {
  
  # resimulate the data
  data <- make_data()
  
  mod <- did::att_gt(yname = "dep_var", 
                     tname = "year",
                     idname = "unit",
                     gname = "cohort_year",
                     control_group = "notyettreated",
                     bstrap = FALSE,
                     data = data,
                     print_details = FALSE)
  event_std <- did::aggte(mod, type = "dynamic")
  
  att.egt <- event_std$att.egt
  names(att.egt) <- event_std$egt
  
  # grab the obs we need
  broom::tidy(att.egt) %>% 
    filter(names %in% -5:5) %>% 
    mutate(t = -5:5, estimate = x) %>% 
    select(t, estimate)
}

data_CS <- map_dfr(1:nrep, run_CS)

ES_plot_CS <- data_CS %>% 
  group_by(t) %>% 
  summarize(avg = mean(estimate),
            sd = sd(estimate),
            lower.ci = avg - 1.96*sd,
            upper.ci = avg + 1.96*sd) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
  ggplot(aes(x = t, y = avg)) + 
  #geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
  geom_point(color = 'blue', size = 3) + 
  geom_line(aes(color = 'Estimated Effect'), size = 1) + 
  geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle("Event-study-parameters estimated using Callaway and Sant'Anna (2021)\nComparison group: Not-yet-treated")+
  scale_color_manual(values = colors) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom", 
        legend.title = element_blank())

ES_plot_CS

#ggsave("es_plot_CS.png", ES_plot_CS, width = 10, height = 5, dpi = 150)