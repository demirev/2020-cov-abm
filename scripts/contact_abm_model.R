source("R/model_abm_functions.R")
source("R/plot_abm_functions.R")

plan("multisession")

# setup -------------------------------------------------------------------
# load contacts data
contacts <- bind_rows(
  read_csv("data/contact_distributions_u18.csv") %>%
    mutate(age = "student"),
  read_csv("data/contact_distributions_o18.csv") %>%
    mutate(age = "adult"),
  read_csv("data/contact_distributions_o18.csv") %>%
    mutate(age = "middle_age"),
  read_csv("data/contact_distributions_o18.csv") %>%
    mutate(age = "pensioner") %>%
    mutate(e_work = 0) # pensioners don't work
)

# define initial population
init_pop <- generate_population(
  household_distribution = list(
    tibble(n = 100000, student = 0, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 100000, student = 1, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 10000, student = 2, adult = 2, middle_age = 2, pensioner = 2),
    tibble(n = 50000, student = 0, adult = 0, middle_age = 0, pensioner = 2)
  ),
  average_workplace_size = 40,
  average_classroom_size = 25
)

# > nrow(init_pop)
# [1] 1100000
# > table(init_pop$age)
# adult middle_age  pensioner    student 
# 620000      50000     190000     240000

# define initial infected individuals
init_inf <- generate_initial_infected(
  init_pop,
  n_initial_infections = c(
    "student" = 5,
    "adult" = 5,
    "middle_age" = 5,
    "pensioner" = 5    
  ),
  contact_distribution = contacts
)

# run simulations ---------------------------------------------------------

# 0. evaluate a baseline policy to test against
simulation_baseline <- simulate_pandemic_policy_sequence(
  initial_population = init_pop,
  initial_infected = init_inf,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = list(list(scenario = "no_measures", n_days = 120)),
  contact_distribution = contacts
)

p0 <- plot_single_simulation(simulation_baseline)

# 1. evaluate all simple strategies and compare them
simulation_all_scenarios <- c(
  "isolation_only","hh_quaratine_only","hh_work_only",
  "isolation_manual_tracing_met_only","isolation_manual_tracing_met_limit",
  "isolation_manual_tracing","cell_phone","cell_phone_met_limit",
  "pop_testing","pt_extra"
) %>%
  future_map(function(scn) {
    simulate_pandemic_policy_sequence(
      initial_population = init_pop,
      initial_infected = init_inf,
      initial_recovered = tibble(individual_id = "")[0,],
      policy_sequence = list(
        list(scenario = "no_measures", n_days = 20), # assume 20 days before measures are taken
        list(scenario = scn, n_days = 100) # let simulation run for total of 120 days
      ),
      contact_distribution = contacts
    )
  })

names(simulation_all_scenarios) <- c(
  "isolation_only","hh_quaratine_only","hh_work_only",
  "isolation_manual_tracing_met_only","isolation_manual_tracing_met_limit",
  "isolation_manual_tracing","cell_phone","cell_phone_met_limit",
  "pop_testing","pt_extra"
)

pA <- plot_different_scenarios(simulation_all_scenarios)

# 3. simulate a complex policy regime many times 
policy <- list(
  list(scenario = "no_measures", n_days = 10),
  list(scenario = "hh_quaratine_only", n_days = 20),
  list(scenario = "hh_work_only", n_days = 20),
  list(scenario = "isolation_manual_tracing_met_only", n_days = 20),
  list(scenario = "cell_phone", n_days = 20),
  list(scenario = "pop_testing", n_days = 30)
)

simulation_specific_policy <- simulate_pandemic_policy_sequence_ntimes(
  initial_population = init_pop,
  initial_infected = init_inf,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = policy,
  n_times = 8, # run 8 simulations of the same policy
  contact_distribution = contacts
)

p1 <- plot_multiple_simulations(simulation_specific_policy)


