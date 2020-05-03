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
    #tibble(n = 50000, student = 2, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 10000, student = 2, adult = 2, middle_age = 2, pensioner = 2),
    #tibble(n = 10000, student = 2, adult = 0, middle_age = 0, pensioner = 2),
    #tibble(n = 50000, student = 0, adult = 0, middle_age = 0, pensioner = 1),
    #tibble(n = 30000, student = 0, adult = 0, middle_age = 1, pensioner = 0),
    #tibble(n = 100000, student = 0, adult = 1, middle_age = 0, pensioner = 0),
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

# define intervention policy to test
policy <- list(
  list(scenario = "no_measures", n_days = 40),
  list(scenario = "hh_quaratine_only", n_days = 20),
  list(scenario = "hh_work_only", n_days = 20),
  list(scenario = "isolation_manual_tracing_met_only", n_days = 20),
  list(scenario = "cell_phone", n_days = 20)
)

# simulate 20 runs of this policy sequence
simulation_results <- simulate_pandemic_policy_sequence_ntimes(
  initial_population = init_pop,
  initial_infected = init_inf,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = policy,
  n_times = 6,
  contact_distribution = contacts
)

# baseline of no measures for comparison - note n_days equals total n_days of target policy
simulation_baseline <- simulate_pandemic_policy_sequence_ntimes(
  initial_population = init_pop,
  initial_infected = init_inf,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = list(list(scenario = "no_measures", n_days = 120)),
  n_times = 1,
  contact_distribution = contacts
)

tictoc::tic()
simulate_pandemic_policy_sequence(
  initial_population = init_pop,
  initial_infected = init_inf,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = list(list(scenario = "no_measures", n_days = 120)),
  contact_distribution = contacts
)
tictoc::toc()

# plot results ------------------------------------------------------------

