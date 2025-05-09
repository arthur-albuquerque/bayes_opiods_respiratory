pacman::p_load(MetaStan, rio, here, dplyr, tidyr, metafor,
               data.table, ggdist, ggplot2, MetBrewer)

# Data ----

dat = import(here("data/raw_data.xlsx"))

dat_OR = escalc(measure="OR",
                ai=Int_Events,        # ai = events in group 1
                n1i=Int_Total,      # n1i = total in group 1
                ci=Control_Events,        # ci = events in group 2
                n2i=Control_Total,      # n2i = total in group 2
                data=dat) |> 
  data.frame()

dat_OR$sei = sqrt(dat_OR$vi)

setnames(dat_OR, old=c("Int_Events","Control_Events"), new=c("r2", "r1"))
setnames(dat_OR, old=c("Int_Total","Control_Total"), new=c("n2", "n1"))

dat_MetaStan <- create_MetaStan_dat(dat = dat_OR,
                                    armVars = c(responders = "r",
                                                sampleSize = "n"))

# Model ----

meta_overall = 
  meta_stan(data = dat_MetaStan,
            seed = 123,
            chains = 4,
            iter = 4000,
            warmup = 2000,
            likelihood = "binomial",
            param = "Smith",               # model parameterization
            re = TRUE,
            ncp = TRUE,
            mu_prior = c(0, 10),           # mu = baseline risk
            theta_prior = c(0, 2.82),      # theta = treatment effect
            tau_prior_dist = "half-normal",
            tau_prior = 0.5,               # tau = heterogeneity
            refresh = 0 # suppress Stan sampling text
  )

meta_subgroup = 
  meta_stan(data = dat_MetaStan,
            seed = 123,
            chains = 4,
            iter = 4000,
            warmup = 2000,
            likelihood = "binomial",
            param = "Smith",               # model parameterization
            re = TRUE,
            ncp = TRUE,
            mu_prior = c(0, 10),           # mu = baseline risk
            theta_prior = c(0, 2.82),      # theta = treatment effect
            tau_prior_dist = "half-normal",
            tau_prior = 0.5,               # tau = heterogeneity
            mreg = TRUE,
            cov = dat_OR$Low_Dose,
            beta_prior = c(0, 10),        # beta = meta-regression
            refresh = 0 # suppress Stan sampling text
            )

# Posterior Samples ----

## Overall ----

log_OR_overall_samples <- rstan::extract(meta_overall$fit)$theta %>% as.vector()

overall_hdi <- ggdist::mean_hdi(log_OR_overall_samples)

overall_summary = 
  list(
    mu_mean = overall_hdi$y,
    mu_lower = overall_hdi$ymin,
    mu_upper = overall_hdi$ymax
  )

# Control risk (mu)
mu_logodds_matrix = rstan::extract(meta_overall$fit)$mu 

# Compute average mu across all studies for each iteration
mu_logodds = rowMeans(mu_logodds_matrix)

model_mean_control_risk = mean(mu_logodds) |> plogis()

## Subgroups ----

log_OR_group_samples <- rstan::extract(meta_subgroup$fit)$theta |>as.vector()
log_OR_beta_samples <- rstan::extract(meta_subgroup$fit)$beta |> as.vector()
log_OR_tau_samples <- rstan::extract(meta_subgroup$fit)$tau |> as.vector()

samples <- tibble(high_dose = log_OR_group_samples,
                  beta = log_OR_beta_samples,
                  low_dose = log_OR_group_samples + log_OR_beta_samples) |> 
  pivot_longer(cols = c(high_dose, low_dose)) |>
  group_by(name)

subgroups_hdi <- samples |> ggdist::mean_hdi(value)
tau_hdi = median_hdi(log_OR_tau_samples)

subgroups_summary = 
  list(
    high_dose_mean = subgroups_hdi$value[subgroups_hdi$name == "high_dose"],
    high_dose_lower = subgroups_hdi$.lower[subgroups_hdi$name == "high_dose"],
    high_dose_upper = subgroups_hdi$.upper[subgroups_hdi$name == "high_dose"],
    
    low_dose_mean = subgroups_hdi$value[subgroups_hdi$name == "low_dose"],
    low_dose_lower = subgroups_hdi$.lower[subgroups_hdi$name == "low_dose"],
    low_dose_upper = subgroups_hdi$.upper[subgroups_hdi$name == "low_dose"],
    
    tau_median = tau_hdi$y,
    tau_lower = tau_hdi$ymin,
    tau_upper = tau_hdi$ymax
  )

## Altogether ----

overall_text = paste0(exp(overall_hdi$y) |> round(2),
                      " [",
                      exp(overall_hdi$ymin) |> round(2),
                      ", ",
                      exp(overall_hdi$ymax) |> round(2),
                      "]")

high_text = paste0(exp(subgroups_summary$high_dose_mean) |> round(2),
                   " [",
                   exp(subgroups_summary$high_dose_lower) |> round(2),
                   ", ",
                   exp(subgroups_summary$high_dose_upper) |> round(2),
                   "]")

low_text = paste0(exp(subgroups_summary$low_dose_mean) |> round(2),
                  " [",
                  exp(subgroups_summary$low_dose_lower) |> round(2),
                  ", ",
                  exp(subgroups_summary$low_dose_upper) |> round(2),
                  "]")

samples = tibble(
  `Overall` = log_OR_overall_samples,
  `High Dose\nSubgroup` = log_OR_group_samples,
  `Low Dose\nSubgroup` = log_OR_group_samples + log_OR_beta_samples
) |> 
  pivot_longer(Overall:`Low Dose\nSubgroup`)

# Posterior Probabilities ----

estimated_risk_fun = function(OR, Rc = model_mean_control_risk){
  # Equation 8 in https://doi.org/10.1016/j.jclinepi.2020.08.019
  (Rc*OR)/(Rc*(OR - 1) + 1)
}

overall_RD_samples = 
 estimated_risk_fun(OR = exp(filter(samples, name == "Overall")$value)) -
  model_mean_control_risk


overall_prob1 = 
  (mean(filter(samples, name == "Overall")$value > log(1)) * 100) |> round(1)

overall_prob2 = 
  (mean(filter(samples, name == "Overall")$value > log(2)) * 100) |> round(1)

low_prob1 = 
  (mean(filter(samples, name == "Low Dose\nSubgroup")$value > log(1)) * 100) |> round(1)

low_prob2 = 
  (mean(filter(samples, name == "Low Dose\nSubgroup")$value > log(2)) * 100) |> round(1)

high_prob1 = 
  (mean(filter(samples, name == "High Dose\nSubgroup")$value > log(1)) * 100) |> round(1)

high_prob2 = 
  (mean(filter(samples, name == "High Dose\nSubgroup")$value > log(2)) * 100) |> round(1)

# Plot ----

# colors
cols = met.brewer(name="VanGogh1", n=7, type="discrete")

# Plot!
samples |> 
  ggplot() +
  aes(x = value, y = name, fill = name,) +
  stat_halfeye(.width = 0.95,
               point_interval = mean_hdi) +
  scale_fill_manual(values = rev(c(cols[2], cols[5], cols[7]))) +
  
  geom_text(
    data = data.frame(
      name = c("Overall", "Low Dose\nSubgroup", "High Dose\nSubgroup"),
      label = c(overall_text, low_text, high_text),  
      x = log(15)  
    ),
    aes(x = x, y = name, label = label),
    hjust = 0,  # Align text to the left (so it extends rightward)
    size = 4,
    nudge_y = 0.15  # Slight vertical nudge to align with distributions
  ) +
  
  geom_text(
    data = data.frame(
      name = c("Overall"),
      label = "Probability of Harm (%)",  # Replace with your desired text
      x = log(120) 
    ),
    aes(x = x, y = name, label = label),
    hjust = 0,  # Align text to the left (so it extends rightward)
    size = 4,
    nudge_y = 0.9  # Slight vertical nudge to align with distributions
  ) +
  
  geom_text(
    data = data.frame(
      name = c("Overall"),
      label = "OR > 1",  # Replace with your desired text
      x = log(200) 
    ),
    aes(x = x, y = name, label = label),
    hjust = 0,  # Align text to the left (so it extends rightward)
    size = 4,
    nudge_y = 0.6  # Slight vertical nudge to align with distributions
  ) +
  
  geom_text(
    data = data.frame(
      name = c("Overall"),
      label = "OR > 2",  # Replace with your desired text
      x = log(900) 
    ),
    aes(x = x, y = name, label = label),
    hjust = 0,  # Align text to the left (so it extends rightward)
    size = 4,
    nudge_y = 0.6  # Slight vertical nudge to align with distributions
  ) +
  
  geom_text(
    data = data.frame(
      name = c("Overall", "Low Dose\nSubgroup", "High Dose\nSubgroup"),
      label = c(overall_prob1, low_prob1, high_prob1),  # Replace with your desired text
      x = log(250)
    ),
    aes(x = x, y = name, label = label),
    hjust = 0,  # Align text to the left (so it extends rightward)
    size = 4,
    nudge_y = 0.15  # Slight vertical nudge to align with distributions
  ) +
  
  geom_text(
    data = data.frame(
      name = c("Overall", "Low Dose\nSubgroup", "High Dose\nSubgroup"),
      label = c(overall_prob2, low_prob2, high_prob2),  # Replace with your desired text
      x = log(1100)
    ),
    aes(x = x, y = name, label = label),
    hjust = 0,  # Align text to the left (so it extends rightward)
    size = 4,
    nudge_y = 0.15  # Slight vertical nudge to align with distributions
  ) +
  
  labs(x = "Odds Ratio (log scale)",
       y = " ") +
  scale_x_continuous(
    breaks = c(1, 2, 5, 10) |> log(),
    labels = function(x) exp(x)
  ) +
  coord_cartesian(x = c(0.6, 2000) |> log()) +
  theme(
    legend.position = "none",
    plot.title.position = 'plot',
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    panel.background = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.3)
  )


# Forest plot ----

forest_data_OR = metabin(
  data = dat_OR,
  studlab = Study,
  
  # Very important: no frequentist analysis in this plot
  subgroup = Low_Dose,
  print.subgroup.name = FALSE,
  common = FALSE,
  random = FALSE,
  
  event.e = r2,  # event.e = events in group 1
  n.e = n2,     # n.e = total in group 1
  event.c = r1,  # event.c = events in group 2
  n.c = n1,     # n.c = total in group 2
  
  sm = "OR"
)

m1 <- metaadd(
  forest_data_OR,
  type = "random", text = "Overall",
  TE = overall_summary$mu_mean,
  lower = overall_summary$mu_lower,
  upper = overall_summary$mu_upper)

# Manually change subgroup results in meta-analysis object
#
mu_mean_subgroup1 <- subgroups_summary$high_dose_mean
mu_lower_subgroup1 <- subgroups_summary$high_dose_lower
mu_upper_subgroup1 <- subgroups_summary$high_dose_upper
#
mu_mean_subgroup2 <- subgroups_summary$low_dose_mean
mu_lower_subgroup2 <- subgroups_summary$low_dose_lower
mu_upper_subgroup2 <- subgroups_summary$low_dose_upper
#
m1$TE.random.w <- c(mu_mean_subgroup1, mu_mean_subgroup2)
m1$lower.random.w <- c(mu_lower_subgroup1, mu_lower_subgroup2)
m1$upper.random.w <- c(mu_upper_subgroup1, mu_upper_subgroup2)
#
names(m1$TE.random.w) <- names(m1$lower.random.w) <-
  names(m1$upper.random.w) <- c("High Dose", "Low Dose")

forest(m1,
       layout = "RevMan5",
       sortvar = TE,
       text.random = "Bayesian Summary",
       smlab = "Odds Ratio (log scale)",
       leftcols=c("studlab", "event.e", "n.e",
                  "event.c", "n.c", "effect.ci"),
       leftlabs=c("Study",  NA, NA, NA, NA, "OR"),
       
       pooled.events = TRUE,
       lab.e = "Experimental", 
       lab.c = "Control",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       colgap = "5mm",
       
       # do not print heterogeneity statistics within subgroups
       hetstat = FALSE,
       # do not print results of test for subgroup differences
       test.subgroup = FALSE)
