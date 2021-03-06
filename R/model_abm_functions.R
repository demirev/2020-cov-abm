library(dplyr)
library(readr)
library(purrr)
library(furrr)
library(tidyr)
library(docstring)

plan(multicore(workers = max(1, availableCores() - 1)))

sample_n_safe <- function(data, n) {
  #' helper - sample_n but check if enough data
  if (nrow(data) < n) return(data) else sample_n(data, n)
}

generate_population <- function(
  household_distribution = list(
    tibble(n = 800, student = 0, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 800, student = 1, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 300, student = 2, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 300, student = 3, adult = 2, middle_age = 0, pensioner = 0),
    tibble(n = 50, student = 1, adult = 2, middle_age = 1, pensioner = 0),
    tibble(n = 50, student = 2, adult = 2, middle_age = 2, pensioner = 0),
    tibble(n = 50, student = 1, adult = 2, middle_age = 1, pensioner = 1),
    tibble(n = 50, student = 2, adult = 2, middle_age = 2, pensioner = 2),
    tibble(n = 50, student = 0, adult = 0, middle_age = 2, pensioner = 0),
    tibble(n = 50, student = 1, adult = 0, middle_age = 2, pensioner = 0),
    tibble(n = 50, student = 2, adult = 0, middle_age = 2, pensioner = 0),
    tibble(n = 50, student = 3, adult = 0, middle_age = 2, pensioner = 0),
    tibble(n = 50, student = 0, adult = 1, middle_age = 2, pensioner = 0),
    tibble(n = 50, student = 2, adult = 0, middle_age = 2, pensioner = 1),
    tibble(n = 50, student = 2, adult = 0, middle_age = 2, pensioner = 2),
    tibble(n = 50, student = 1, adult = 0, middle_age = 0, pensioner = 2),
    tibble(n = 50, student = 2, adult = 0, middle_age = 0, pensioner = 2),
    tibble(n = 50, student = 2, adult = 0, middle_age = 0, pensioner = 2),
    tibble(n = 50, student = 0, adult = 0, middle_age = 0, pensioner = 1),
    tibble(n = 50, student = 0, adult = 0, middle_age = 1, pensioner = 0),
    tibble(n = 800, student = 0, adult = 1, middle_age = 0, pensioner = 0)
  ),
  average_workplace_size = 30,
  average_classroom_size = 20
) {
  #' generate a tibble with all agents
  #' @example {init_pop <- generate_population()}
  
  population <- household_distribution %>%
    seq_along() %>%
    future_map(function(i) {
      dist <- household_distribution[[i]]
      tibble(
        household_id = rep(
          paste0(i, "_", 1:dist$n), 
          each = dist$student + dist$adult + dist$middle_age + dist$pensioner
        )
      ) %>%
        mutate(
          individual_id = paste0(
            i, "_", household_id, "_",
            1:nrow(.)
          ),
          age = rep(
            c(
              rep("student", dist$student), 
              rep("adult", dist$adult),
              rep("middle_age", dist$middle_age),
              rep("pensioner", dist$pensioner)
            ),
            dist$n
          )
        )
    }) %>%
    reduce(bind_rows) %>%
    mutate(
      workplace_id = ifelse(
        !age %in% c("adult", "middle_age"),
        NA,
        sample(
          1:round(nrow(filter(., age %in% c("adult", "middle_age")))/average_workplace_size), 
          size = nrow(filter(., age %in% c("adult", "middle_age"))),
          replace = T
        )
      ), # id for workplace
      workplace_id = ifelse(
        age %in% c("adult", "middle_age"),
        workplace_id,
        sample(
          1:round(nrow(filter(., !age %in% c("adult", "middle_age")))/average_classroom_size), 
          size = nrow(filter(., !age %in% c("adult", "middle_age"))),
          replace = T
        )
      ), # id for school
      workplace_id = ifelse(
        age %in% c("adult", "middle_age"),
        paste0("W_", workplace_id),
        paste0("S_", workplace_id) # differentiate between both
      )
    )
}

generate_initial_infected <- function(
  population,
  n_initial_infections = c(
    "student" = 3,
    "adult" = 3,
    "middle_age" = 3,
    "pensioner" = 3    
  ),
  contact_distribution = bind_rows(
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
) {
  #' sample initial infected population
  #' @example {init_inf <- generate_initial_infected(init_pop)}
  
  n_initial_infections %>%
    map2(names(.), function(n_i, a) {
      population %>%
        filter(age == a) %>%
        sample_n(n_i) %>%
        bind_cols(
          contact_distribution %>%
            filter(age == a) %>%
            sample_n(n_i)
        )
    }) %>%
    reduce(bind_rows) %>%
    mutate(infected_since = 0)
}

get_all_parameters <- function(
  scenario,
  inf_period,
  p_pop_test,
  cc_risk,
  hh_risk,
  wfh_prob,
  pt_extra,
  pt_extra_reduce,
  trace_adherence,
  isolate_distn,
  prob_symp,
  prob_t_asymp,
  trace_hh,
  trace_prop,
  max_low_fix,
  app_cov,
  met_before_w,
  met_before_h,
  met_before_o
) {
  #' helper to make the following two functions 'purer'
  
  params <- list()
  
  params$inf_period <- inf_period
  params$p_pop_test <- p_pop_test
  params$cc_risk <- cc_risk
  params$hh_risk <- hh_risk
  params$wfh_prob <- wfh_prob
  params$pt_extra <- pt_extra
  params$pt_extra_reduce <- pt_extra_reduce
  
  # Symptomatic and proportion getting tested
  params$trace_adherence <- trace_adherence
  params$p_tested <- trace_adherence # Proportion who get tested
  params$time_isolate <- isolate_distn # Distribution over symptomatic period
  params$p_symptomatic <- prob_symp
  params$transmission_asymp <- prob_t_asymp
  params$phone_coverage <- 1 # App coverage in non-app scenarios
  
  # Set contact limit default high to avoid censoring in default scenarios
  params$max_contacts <- 2e3 
  
  # Tracing parameters - need to reset
  params$hh_trace <- trace_hh # Tracing in HH
  params$ww_trace <- trace_prop # Tracing at work
  params$other_trace <- trace_prop # Tracing others
  
  params$met_before_w <- met_before_w
  params$met_before_h <- met_before_h
  params$met_before_o <- met_before_o
  
  # Define scenario parameters
  if (scenario == "no_measures") {
    params$do_isolation <- F
    params$do_tracing <- F
  }
  
  if (scenario == "isolation_only") {
    params$do_isolation <- T
    params$do_tracing <- F
  }
  
  if (scenario == "hh_quaratine_only") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$ww_trace <- 0 # Tracing at work
    params$other_trace <- 0 # Tracing others
  }
  
  if (scenario == "hh_work_only") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$met_before_w <- 1 # Met before
    params$met_before_o <- 1 # Met before
    params$other_trace <- 0 # Tracing others
  }
  
  if (scenario == "isolation_manual_tracing_met_limit") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$max_contacts <- max_low_fix
  }
  
  if (scenario == "isolation_manual_tracing_met_only") {
    params$do_isolation <- T
    params$do_tracing <- T
  }
  
  if (scenario == "isolation_manual_tracing") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$met_before_w <- 1 # Met before
    params$met_before_o <- 1 # Met before
  }
  
  if (scenario == "cell_phone") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$met_before_w <- 1 # Met before
    params$met_before_o <- 1 # Met before
    params$phone_coverage <- app_cov
  }
  
  if (scenario == "cell_phone_met_limit") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$met_before_w <- 1 # Met before
    params$met_before_o <- 1 # Met before
    params$phone_coverage <- app_cov
    params$max_contacts <- max_low_fix
  }
  
  if (scenario == "pop_testing") {
    params$do_isolation <- F
    params$do_tracing <- F
  }
  
  if (scenario == "pt_extra") {
    params$do_isolation <- T
    params$do_tracing <- T
    params$ww_trace <- 0 # Tracing at work
    params$other_trace <- 0 # Tracing others
  }
  
  params
}

generate_per_person_parameters <- function(
  infected,
  params,
  scenario = "no_measures"
) {
  #' all parameters for the infected population
  
  infected %>% 
    bind_rows(
      tibble(
        wfh = F, phone = F, symp = F, tested = F, inf_period = 5,
        met_before_w = 1, extra_red = 0
      )[0,] # make sure all columns are present
    ) %>%
    mutate( # for new cases (infected_since == 0) draw parameter values
      wfh = infected_since == 0 & runif(nrow(infected)) < params$wfh_prob,
      wfh = ifelse(age == "student", F, wfh), # TODO introdue school vacations
      phone = infected_since == 0 & runif(nrow(infected)) < params$phone_coverage,
      symp = infected_since == 0 & runif(nrow(infected)) < params$p_symptomatic,
      tested = infected_since == 0 & runif(nrow(infected)) < params$p_tested,
      inf_period = ifelse(
        infected_since == 0 & symp & params$do_isolation,
        sample(0:params$inf_period, nrow(infected), prob = params$time_isolate, replace = T),
        ifelse(
          infected_since == 0,
          params$inf_period,
          inf_period
        )
      ),
      met_before_w = ifelse(age == "student", 0.9, params$met_before_w)
    ) %>%
    mutate( # additional adjustments if pop_testing
      pop_tested_ = (scenario == "pop_testing" & 
        infected_since == 0 & 
        runif(nrow(infected)) < params$p_pop_test),
      tested = ifelse(
        pop_tested_,
        T,
        tested
      ),
      inf_period = ifelse(
        pop_tested_,
        sample(c(0:params$inf_period,5), nrow(infected), replace = T),
        inf_period
      )
    ) %>%
    select(-pop_tested_) %>%
    mutate( # additional adjustments if pt_extra
      pt_extra_ = (scenario=="pt_extra" & 
        runif(nrow(infected)) < params$pt_extra),
      tested = ifelse(pt_extra_, T, tested),
      extra_red = ifelse(pt_extra_,  (1 - params$pt_extra_reduce), 1)
    ) %>%
    select(-pt_extra_) %>%
    mutate( # if not covered above then false
      symp = ifelse(is.na(symp), F, symp), 
      tested = ifelse(is.na(tested), F, tested),
      wfh = ifelse(is.na(wfh), F, wfh)
    )
}

simulate_per_person_infections <- function( 
  infected,
  params,
  population = population,
  recovered = recovered,
  scenario = "no_measures",
  household_sizes = population %>%
    group_by(household_id) %>%
    count() %>%
    mutate(home_c = n - 1) # don't count the individual doing the disease spreading
) {
  #' all contacts and infections due to a cohort of infected individuals
  
  result <- infected %>%
    mutate( # Set relative transmission of asymptomatics
      inf_propn = ifelse(symp, 1, params$transmission_asymp) 
    ) %>%
    mutate( # Check if contacts phone traced (in cell phne scenario):
      ww_trace = ifelse(
        scenario == "cell_phone" | scenario == "cell_phone_met_limit",
        ifelse(phone, params$phone_coverage, 0),
        params$ww_trace
      ),
      other_trace = ifelse(
        scenario == "cell_phone" | scenario == "cell_phone_met_limit",
        ifelse(phone, params$phone_coverage, 0),
        params$other_trace
      )
    ) %>%
    mutate(# Proportion infectious
      inf_ratio = as.numeric(infected_since <= inf_period)
    ) %>%
    left_join(
      select(household_sizes, household_id, home_c),
      by = "household_id"
    ) %>%
    mutate( # Tally contacts
      home_c = ifelse(infected_since > 0, 0, home_c),
      work_c = ifelse(is.na(e_work), 0, e_work),
      other_c = ifelse(is.na(e_other), 0, e_other),
      scale_other = min(1, (params$max_contacts/other_c)) # scale down based on max other contacts
    )

  # Draw ids of contacts
  home_c_all_ids <- population %>%
    right_join(
      result %>% 
        select(household_id, home_c, infected_by = individual_id), 
      by = "household_id"
    ) %>% # all matching houyseholds to infected
    anti_join(
      select(result, household_id, individual_id), 
      by = c("household_id", "individual_id")
    ) %>% # remove the individuals doing the infection themselves
    group_by(household_id) %>%
    #mutate(home_c = ifelse(home_c > n(), n(), home_c)) %>% # make sure we are not sampling more than the available number
    sample_n(home_c, replace = F) %>% # sample contacts
    ungroup() %>%
    anti_join(
      select(infected, individual_id), by = "individual_id"
    ) %>% # remove otherwise infected (! only after sampling contacts)
    anti_join(
      select(recovered, individual_id), by = "individual_id"
    ) %>%
    select(infected_by, individual_id) %>%
    nest(individual_id) %>%
    rename(home_c_id = data, individual_id = infected_by) %>%
    mutate(home_c_id = map(home_c_id, ~.$individual_id))
  
  work_c_all_ids <- population %>%
    right_join(
      result %>% 
        select(workplace_id, work_c, infected_by = individual_id), 
      by = "workplace_id"
    ) %>% # all matching workplaces to infected
    anti_join(
      select(result, workplace_id, individual_id), 
      by = c("workplace_id", "individual_id")
    ) %>% # remove the individuals doing the infection themselves
    group_by(workplace_id) %>%
    mutate(work_c = ifelse(work_c > n(), n(), work_c)) %>% # make sure we are not sampling more than the available number
    sample_n(work_c, replace = F) %>% # sample contacts
    ungroup() %>%
    anti_join(
      select(infected, individual_id), by = "individual_id"
    ) %>% # remove otherwise infected (! only after sampling contacts)
    anti_join(
      select(recovered, individual_id), by = "individual_id"
    ) %>%
    select(infected_by, individual_id) %>%
    nest(individual_id) %>%
    rename(work_c_id = data, individual_id = infected_by) %>%
    mutate(work_c_id = map(work_c_id, ~.$individual_id))
  
  other_c_all_ids <- population %>%
    sample_n(
      sum(result$other_c), 
      replace = T
    ) %>% # sample other encounters from entire population (Note: this allows for meeting yourself and potentially meeting the same person twice but the chances are small and checking against that is too comp. expensive)
    mutate(
      infected_by = rep(result$individual_id, result$other_c)
    ) %>%
    anti_join(
      select(infected, individual_id), by = "individual_id"
    ) %>% # remove otherwise infected (! only after sampling contacts)
    anti_join(
      select(recovered, individual_id), by = "individual_id"
    ) %>%
    select(infected_by, individual_id) %>%
    nest(individual_id) %>%
    rename(other_c_id = data, individual_id = infected_by) %>%
    mutate(other_c_id = map(other_c_id, ~.$individual_id))
  
  
  result <- result %>%
    left_join(home_c_all_ids, by = "individual_id") %>%
    left_join(work_c_all_ids, by = "individual_id") %>%
    left_join(other_c_all_ids, by = "individual_id") %>%
    mutate( # make sure no null values for id var
      home_c_id = ifelse(
        map_lgl(home_c_id, is.null),
        map(nrow(.), ~c("")[0]),
        home_c_id
      ),
      work_c_id = ifelse(
        map_lgl(work_c_id, is.null),
        map(nrow(.), ~c("")[0]),
        work_c_id
      ),
      other_c_id = ifelse(
        map_lgl(other_c_id, is.null),
        map(nrow(.), ~c("")[0]),
        other_c_id
      )
    ) %>%
    mutate( # Generate basic infections
      home_inf_basic = rbinom(
        nrow(infected), 
        map_dbl(home_c_id, length), 
        prob = params$hh_risk * inf_propn
      ),
      work_inf_basic = rbinom(
        nrow(infected), 
        map_dbl(work_c_id, length), 
        prob = params$cc_risk * inf_propn
      ),
      other_inf_basic = rbinom(
        nrow(infected), 
        map_dbl(other_c_id, length), 
        prob = params$cc_risk * inf_propn
      ),
      rr_basic = home_inf_basic + work_inf_basic + other_inf_basic
    ) %>%
    mutate( # Gnerate infections
      inf_ratio_w = ifelse(wfh, 0, inf_ratio), # check if working from home
      home_infect = rbinom(
        nrow(infected), 
        home_inf_basic, 
        prob = inf_ratio
      ),
      work_infect = rbinom(
        nrow(infected), 
        work_inf_basic, 
        prob = inf_ratio_w * extra_red
      ),
      other_infect = rbinom(
        nrow(infected),
        other_inf_basic, 
        prob = inf_ratio * scale_other * extra_red # scale by maximum
      ), 
      rr = home_infect + work_infect + other_infect
    ) %>%
    mutate( # Cntact tracing - tally contacts traced
      home_traced = rbinom(
        nrow(infected), 
        map_dbl(home_c_id, length), 
        prob = params$hh_trace
      ),
      work_traced = rbinom(
        nrow(infected), 
        map_dbl(work_c_id, length), 
        prob = ww_trace * met_before_w
      ), # !met_before_w defined in data not in params
      other_traced = rbinom(
        nrow(infected), 
        map_dbl(other_c_id, length), 
        prob = params$other_trace * params$met_before_o * scale_other # NOTE: different from AK original code but I suspect it was a mistake in the original
      )
    ) %>%
    mutate( # Infections averted
      home_averted = ifelse(
        tested & symp & params$do_tracing,
        rbinom(
          nrow(infected),
          home_infect, # TODO not sure if this should be home_traced instead...
          prob = params$hh_trace * params$trace_adherence # ... and this just trace_adherence?
        ),
        0        
      ),
      work_averted = ifelse(
        tested & symp & params$do_tracing,
        rbinom(
          nrow(infected),
          work_infect,
          prob = ww_trace * met_before_w * params$trace_adherence
        ),
        0
      ),
      other_averted = ifelse(
        tested & symp & params$do_tracing,
        rbinom(
          nrow(infected),
          other_infect,
          prob = params$met_before_o * params$other_trace * params$trace_adherence
        ),
        0
      ),
      total_averted = home_averted + work_averted + other_averted,
      rr_reduced = rr - total_averted
    ) %>%
    mutate( # Final ids of all infected
      home_i_id = pmap(
        list(home_c_id, home_infect, home_averted),
        function(ids, i, a) sample(ids, size = i - a)
      ),
      work_i_id = pmap(
        list(work_c_id, work_infect, work_averted),
        function(ids, i, a) sample(ids, size = i - a)
      ),
      other_i_id = pmap(
        list(other_c_id, other_infect, other_averted),
        function(ids, i, a) sample(ids, size = i - a)
      )
    ) %>%
    mutate( # Count total traced
      total_traced = 0,
      total_traced = ifelse(
        scenario == "hh_quaratine_only",
        home_traced,
        total_traced
      ),
      total_traced = ifelse(
        scenario == "hh_work_only",
        home_traced + work_traced,
        total_traced
      ),
      total_traced = ifelse(
        scenario %in% c(
          "isolation_manual_tracing_met_only", 
          "isolation_manual_tracing_met_limit",
          "isolation_manual_tracing", 
          "cell_phone",
          "cell_phone_met_limit"
        ),
        home_traced + work_traced + other_traced ,
        total_traced
      )
    ) %>%
    mutate( # Format outputs
      home_infected_reduced = home_infect - home_averted,
      work_infected_reduced = work_infect - work_averted,
      other_infected_reduced = other_infect - other_averted,
      infected_id = pmap(
        list(home_i_id, work_i_id, other_i_id),
        function(h,w,o) c(h,w,o)
      ),
      infected_at = pmap(
        list(home_i_id, work_i_id, other_i_id),
        function(h,w,o) {
          c(
            rep("home", length(h)),
            rep("work", length(w)),
            rep("other", length(o))
          )
        }
      )
    ) %>%
    select( # Outputs
      individual_id,
      rr = rr_basic, 
      rr_reduced,
      home_infected_reduced,
      work_infected_reduced,
      other_infected_reduced,
      total_traced = total_traced,
      infected_id,
      infected_at
    ) %>%
    unnest() %>% # unlist infected_id and infected_at columns
    distinct(infected_id, .keep_all = T) # don't double count people infected by more than 1 person
}

infect_daily <- function(
  population,
  infected,
  recovered = tibble(id = "")[0,],
  scenario = "no_measures",
  household_sizes = population %>%
    group_by(household_id) %>%
    count() %>%
    mutate(home_c = n - 1),
  max_low_fix = 4, # Social distancing limit in these scenarios
  wfh_prob = 0, # Probability people are working from home
  trace_prop = 0.95, # Proportion of contacts traced
  app_cov = 0.53, # App coverage
  prob_symp = 0.6, # Proportion symptomatic
  prob_t_asymp = 0.5, # Transmission probability if asymptotic
  isolate_distn = c(0,0.25,0.25,0.2,0.3,0), # distribution of time to isolate (1st day presymp)
  pt_extra = 0, # Optional extra transmission intervention
  pt_extra_reduce = 0, # Reduction from extra intervention
  hh_risk = 0.2, # HH risk
  cc_risk = 0.06, # Outside HH contact risk
  trace_adherence = 0.9, 
  p_pop_test = 0.05, # Proportion mass tested (5% per week)
  inf_period = 5, # Infectious period
  trace_hh = 1,
  met_before_w = 0.79, # At work. At school = 90%, which is defined in function later on
  met_before_h = 1, # Within HH
  met_before_o =  0.52 # In other settings
) {
  #' simulate all contacts for all inidividuals in one day
  #' @example {daily_result <- infect_daily(init_pop, init_inf)}
  
  scenario_list <- c(
    "no_measures","isolation_only","hh_quaratine_only","hh_work_only",
    "isolation_manual_tracing_met_only","isolation_manual_tracing_met_limit",
    "isolation_manual_tracing","cell_phone","cell_phone_met_limit",
    "pop_testing","pt_extra"
  )
  
  if (!scenario %in% scenario_list) {
    stop("invalid scenario")
  }
  
  # bundle all parameters in list so to be easily passable to functions below
  params <- get_all_parameters(
    scenario,
    inf_period = inf_period,
    p_pop_test = p_pop_test,
    cc_risk = cc_risk,
    hh_risk = hh_risk,
    wfh_prob = wfh_prob,
    pt_extra = pt_extra,
    pt_extra_reduce = pt_extra_reduce,
    trace_adherence = trace_adherence,
    isolate_distn = isolate_distn,
    prob_symp = prob_symp,
    prob_t_asymp = prob_t_asymp,
    trace_hh = trace_hh,
    trace_prop = trace_prop,
    max_low_fix = max_low_fix,
    app_cov = app_cov,
    met_before_w = met_before_w,
    met_before_h = met_before_h,
    met_before_o = met_before_o
  ) 
  
  # update infection paramteres
  infected <- infected %>%
    generate_per_person_parameters(params = params, scenario = scenario)
  
  # simulate infections
  daily_r <- infected %>%
    simulate_per_person_infections(
      params = params,
      population = population,
      recovered = recovered,
      scenario = scenario,
      household_sizes = household_sizes
    )
  
  return(list(daily_r = daily_r, infected = infected))
  
}

draw_contact_rates <- function(
  tb,
  contact_distribution = bind_rows(
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
) {
  #' helper - draw contact rates for newly infected
  tb_by_age <- count(tb, age)
  
  if (nrow(tb_by_age)) {
    res <- count(tb, age) %>%
      split(1:nrow(.)) %>%
      map(function(age_group) {
        contact_distribution %>%
          filter(age == age_group$age) %>%
          sample_n(age_group$n, replace = T) %>%
          select(-age) %>%
          bind_cols(filter(tb, age == age_group$age))
      }) %>%
      reduce(bind_rows) 
  } else {
    res <- bind_cols(
      tb_by_age,
      tibble(e_home = 0, e_work = 0, e_other = 0)[0,] # just to have the cols
    )
  }
  
  res
}

simulate_pandemic_days <- function(
  initial_population,
  initial_infected,
  initial_recovered = tibble(individual_id = "")[0,],
  n_days = 100,
  is_initial = T,
  scenario = "no_measures",
  inf_period = 5, # Infectious period
  contact_distribution = bind_rows(
    read_csv("data/contact_distributions_u18.csv") %>%
      mutate(age = "student"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "adult"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "middle_age"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "pensioner") %>%
      mutate(e_work = 0) # pensioners don't work
  ),
  verbose = T,
  # max_low_fix = 4, # Social distancing limit in these scenarios
  # wfh_prob = 0, # Probability people are working from home
  # trace_prop = 0.95, # Proportion of contacts traced
  # n_run = 5e3, # Number of simualtions
  # app_cov = 0.53, # App coverage
  # prob_symp = 0.6, # Proportion symptomatic
  # prob_t_asymp = 0.5, # Transmission probability if asymptotic
  # isolate_distn = c(0,0.25,0.25,0.2,0.3,0), # distribution of time to isolate (1st day presymp)
  # pt_extra = 0, # Optional extra transmission intervention
  # pt_extra_reduce = 0, # Reduction from extra intervention
  # hh_risk = 0.2, # HH risk
  # cc_risk = 0.06, # Outside HH contact risk
  # trace_adherence = 0.9, 
  # p_pop_test = 0.05, # Proportion mass tested (5% per week)
  # trace_hh = 1,
  # met_before_w = 0.79, # At work. At school = 90%, which is defined in function later on
  # met_before_h = 1, # Within HH
  # met_before_o =  0.52 # In other settings
  ...
) {
  #' simulate a number of days since the pandemic outbreak
  #' @example {result <- simulate_pandemic_days(init_pop, init_inf, n_days = 10)}
  
  t_population <- initial_population
  t_infected <- initial_infected
  t_recovered <- initial_recovered
  
  if (is_initial) {
    ts_infected <- tibble(
      individual_id = "",
      infected_id = initial_infected$individual_id,
      infected_at = "",
      infected_since = 0,
      on_day = 0
    )
  } else {
    ts_infected <- tibble()
  } #TODO rembember why I did this... since its immedeately overwritten below
  ts_infected <- tibble()
  ts_recovered <- tibble()
  ts_susceptible <- tibble()
  
  household_sizes <- t_population %>%
    group_by(household_id) %>%
    count() %>%
    mutate(home_c = n - 1) # this doesn't change so calculate once to save time
  
  for (day in seq(n_days)) {
    # get new infections
    daily_infected <- infect_daily(
      population = t_population,
      infected = t_infected,
      recovered = t_recovered,
      scenario = scenario,
      inf_period = inf_period,
      household_sizes = household_sizes,
      ...
    )
    
    # update infected info (some of it gets determined inside infect_daily)
    t_infected <- daily_infected$infected
    
    # age current infections with one day
    t_infected$infected_since <- t_infected$infected_since + 1
    
    # recover some patients
    t_recovered <- bind_rows(
      t_recovered,
      t_infected %>%
          filter(infected_since > inf_period)
    )
    
    # add new infections
    t_infected <- t_infected %>%
      filter(infected_since <= inf_period) 
    
    if (nrow(daily_infected$daily_r)) {
      t_infected <- t_infected %>%
        bind_rows(
          tibble(
            individual_id = daily_infected$daily_r$infected_id,
            infected_since = 0
          ) %>%
            left_join(
              select(
                t_population, individual_id, household_id, workplace_id, age
              ), 
              by = c("individual_id")
            ) %>%
            draw_contact_rates(contact_distribution)
        )
    }
    
    # record flows
    ts_infected <- bind_rows(
      ts_infected, 
      mutate(daily_infected$daily_r, on_day = day)
    )
    ts_recovered <- bind_rows(
      ts_recovered, 
      tibble(n = nrow(t_recovered), on_day = day)
    )
    ts_susceptible <- bind_rows(
      ts_susceptible,
      tibble(n = nrow(t_population) - nrow(t_recovered) - nrow(t_infected), on_day = day)
    )
    
    if (verbose) {
      cat(paste0(nrow(daily_infected$daily_r), " infected on day ", day, "\n"))
    }
  }
  
  return(list(
    "final_population" = t_population,
    "final_infected" = t_infected,
    "final_recovered" = t_recovered,
    "s_ts" = ts_susceptible,
    "i_ts" = ts_infected,
    "r_ts" = ts_recovered
  ))
} 

simulate_pandemic_policy_sequence <- function(
  initial_population,
  initial_infected,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = list(
    list(scenario = "no_measures", n_days = 20)
  ),
  inf_period = 5, # Infectious period
  contact_distribution = bind_rows(
    read_csv("data/contact_distributions_u18.csv") %>%
      mutate(age = "student"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "adult"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "middle_age"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "pensioner") %>%
      mutate(e_work = 0) # pensioners don't work
  ),
  ...
) {
  #' simulate a sequence of different intervention policies
  #' @example {result_sequence <- simulate_pandemic_policy_sequence(
  #' init_pop, init_inf, 
  #' policy_sequence = list(
  #' list(scenario = "no_measures", n_days = 5),
  #' list(scenario = "no_measures", n_days = 5)
  #' )
  #' )
  #' }
  
  t_population <- initial_population
  t_infected <- initial_infected
  t_recovered <- initial_recovered
  
  ts_infected <- tibble(
    individual_id = "",
    infected_id = initial_infected$individual_id,
    infected_at = "",
    on_day = 0
  )
  ts_recovered <- tibble()
  ts_susceptible <- tibble()
  
  days_to_add <- 0
  
  for (policy in policy_sequence) {
    
    # simulat policy for policy duration days
    policy_results <- simulate_pandemic_days(
      initial_population = t_population,
      initial_infected = t_infected,
      initial_recovered = t_recovered,
      n_days = policy$n_days,
      scenario = policy$scenario,
      inf_period = inf_period,
      contact_distribution = contact_distribution, 
      is_initial = F,
      ...
    )
    
    # update SIR info
    t_infected <- policy_results$final_infected
    t_recovered <- policy_results$final_recovered
    
    # record flows
    ts_infected <- bind_rows(
      ts_infected, 
      mutate(policy_results$i_ts, on_day = on_day + days_to_add)
    )
    ts_recovered <- bind_rows(
      ts_recovered, 
      mutate(policy_results$r_ts, on_day = on_day + days_to_add)
    )
    ts_susceptible <- bind_rows(
      ts_susceptible,
      mutate(policy_results$s_ts, on_day = on_day + days_to_add)
    )
    
    days_to_add <- days_to_add + policy$n_days
  }
  
  return(list(
    "final_population" = t_population,
    "final_infected" = t_infected,
    "final_recovered" = t_recovered,
    "s_ts" = ts_susceptible,
    "i_ts" = ts_infected,
    "r_ts" = ts_recovered
  ))
}

simulate_pandemic_policy_sequence_ntimes <- function(
  initial_population,
  initial_infected,
  initial_recovered = tibble(individual_id = "")[0,],
  policy_sequence = list(
    list(scenario = "no_measures", n_days = 20)
  ),
  n_times = 10,
  inf_period = 5, # Infectious period
  contact_distribution = bind_rows(
    read_csv("data/contact_distributions_u18.csv") %>%
      mutate(age = "student"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "adult"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "middle_age"),
    read_csv("data/contact_distributions_o18.csv") %>%
      mutate(age = "pensioner") %>%
      mutate(e_work = 0) # pensioners don't work
  ),
  ...
) {
  #' perform n_times simulations using the same parameters
  #' @example {result_sequence_ntimes <- simulate_pandemic_policy_sequence_ntimes(
  #' init_pop, init_inf, n_times = 3,
  #' policy_sequence = list(
  #' list(scenario = "no_measures", n_days = 5),
  #' list(scenario = "no_measures", n_days = 5)
  #' )
  #' )
  #' }
  
  result <- seq(n_times) %>%
    future_map(
      simulate_pandemic_policy_sequence,
      initial_population = initial_population,
      initial_infected = initial_infected,
      initial_recovered = initial_recovered,
      policy_sequence = policy_sequence,
      inf_period = inf_period,
      contact_distribution = contact_distribution,
      ...
    )
  names(result) = paste0("simulation", seq(n_times))
  result
}
