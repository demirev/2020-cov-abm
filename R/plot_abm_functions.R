library(ggplot2)
library(dplyr)
library(docstring)

plot_single_simulation <- function(simulation_result) {
  #' all plots for a single simulation run
  sir_data <- simulation_result$i_ts %>%
    select(infected_id, on_day) %>%
    distinct(infected_id, .keep_all = T) %>%
    count(on_day) %>%
    mutate(group = "I") %>%
    bind_rows(mutate(simulation_result$s_ts, group = "S")) %>%
    bind_rows(mutate(simulation_result$r_ts, group = "R"))
  
  i_plot <- sir_data %>%
    filter(group == "I") %>%
    ggplot(aes(x = on_day, y = n, color = group)) +
    geom_line() +
    theme_bw() +
    xlab('day') +
    ggtitle("Number Infected By Day")
  
  sir_plot <- sir_data %>%
    ggplot(aes(x = on_day, y = n, color = group)) +
    geom_line() +
    theme_bw() +
    xlab('day') +
    ggtitle("SIR By Day")
  
  rt_data <- simulation_result$i_ts %>%
    select(infected_id, on_day) %>%
    distinct(infected_id, .keep_all = T) %>%
    rename(individual_id = infected_id) %>%
    left_join(
      simulation_result$i_ts %>%
        filter(individual_id != "") %>%
        select(individual_id, infected_id) %>%
        distinct() %>%
        count(individual_id),
      by = "individual_id"
    ) %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  # fivenum(rt_data$n)
  
  rt_plot <- rt_data %>%
    group_by(on_day) %>%
    summarize(n = mean(n)) %>%
    ggplot(aes(x = on_day, y = n)) +
    geom_line() +
    theme_bw() +
    xlab('day') +
    ggtitle("R_t by day") +
    geom_hline(yintercept = 1, linetype = 'dashed')
  
  list(
    sir_data = sir_data,
    rt_data = rt_data,
    sir_plot = sir_plot,
    i_plot = i_plot,
    rt_plot = rt_plot
  )
}

plot_different_scenarios <- function(simulation_results) {
  #' plot and compare different scenarios
  
  sir_data <- map2(simulation_results, names(simulation_results), function(simulation_result, nm) {
    simulation_result$i_ts %>%
      select(infected_id, on_day) %>%
      distinct(infected_id, .keep_all = T) %>%
      count(on_day) %>%
      mutate(group = "I") %>%
      bind_rows(mutate(simulation_result$s_ts, group = "S")) %>%
      bind_rows(mutate(simulation_result$r_ts, group = "R")) %>%
      mutate(scenario = nm)
  }) %>%
    reduce(bind_rows)
  
  i_plot <- sir_data %>%
    filter(group == "I") %>%
    ggplot(aes(x = on_day, y = n, color = scenario, linetype = scenario)) +
    geom_line() +
    theme_bw() +
    xlab('day') +
    ggtitle("Number Infected By Day")
  
  rt_data <- map2(simulation_results, names(simulation_results), function(simulation_result, nm) {
    simulation_result$i_ts %>%
      select(infected_id, on_day) %>%
      distinct(infected_id, .keep_all = T) %>%
      rename(individual_id = infected_id) %>%
      left_join(
        simulation_result$i_ts %>%
          filter(individual_id != "") %>%
          select(individual_id, infected_id) %>%
          distinct() %>%
          count(individual_id),
        by = "individual_id"
      ) %>%
      mutate(n = ifelse(is.na(n), 0, n)) %>%
      mutate(scenario = nm)
  }) %>%
    reduce(bind_rows)
  
  # fivenum(rt_data$n)
  
  rt_plot <- rt_data %>%
    group_by(on_day, scenario) %>%
    summarize(n = mean(n)) %>%
    ungroup() %>%
    filter(on_day <= 0.95*max(on_day)) %>%
    ggplot(aes(x = on_day, y = n, linetype = scenario, color = scenario)) +
    geom_line() +
    theme_bw() +
    xlab('day') +
    ggtitle("R_t by day") +
    geom_hline(yintercept = 1)
  
  list(
    sir_data = sir_data,
    rt_data = rt_data,
    i_plot = i_plot,
    rt_plot = rt_plot
  )
}

plot_multiple_simulations <- function(simulation_results) {
  #' average results of multiple simulatiuon runs
  sir_data <- imap(simulation_results, function(simulation_result, ind) {
    simulation_result$i_ts %>%
      select(infected_id, on_day) %>%
      distinct(infected_id, .keep_all = T) %>%
      count(on_day) %>%
      mutate(group = "I") %>%
      bind_rows(mutate(simulation_result$s_ts, group = "S")) %>%
      bind_rows(mutate(simulation_result$r_ts, group = "R")) %>%
      mutate(simulation = ind)
  }) %>%
    reduce(bind_rows)
  
  sir_data <- bind_rows(
    sir_data,
    sir_data %>%
      group_by(on_day, group) %>%
      summarise(n = mean(n)) %>%
      ungroup() %>%
      mutate(simulation = 'average')
  )
  
  i_plot <- sir_data %>%
    filter(group == "I") %>%
    mutate(type = ifelse(simulation == "average", "average", "simulations")) %>%
    ggplot(aes(x = on_day, y = n, group = simulation, color = type, alpha = type)) +
    geom_line() +
    scale_alpha_manual(values = c(1, 0.5)) +
    scale_color_manual(values = c('black', 'gray80')) +
    theme_bw() +
    xlab('day') +
    ggtitle("Number Infected By Day")
  
  sir_plot <- sir_data %>%
    filter(simulation == "average") %>%
    ggplot(aes(x = on_day, y = n, color = group)) +
    geom_line() +
    theme_bw() +
    xlab('day') +
    ggtitle("SIR By Day")
  
  rt_data <- imap(simulation_results, function(simulation_result, ind) {
    simulation_result$i_ts %>%
      select(infected_id, on_day) %>%
      distinct(infected_id, .keep_all = T) %>%
      rename(individual_id = infected_id) %>%
      left_join(
        simulation_result$i_ts %>%
          filter(individual_id != "") %>%
          select(individual_id, infected_id) %>%
          distinct() %>%
          count(individual_id),
        by = "individual_id"
      ) %>%
      mutate(n = ifelse(is.na(n), 0, n)) %>%
      mutate(simulation = ind)
  }) %>%
    reduce(bind_rows)
  
  # fivenum(rt_data$n)
  
  rt_plot <- rt_data %>%
    group_by(on_day, simulation) %>%
    summarize(n = mean(n)) %>%
    bind_rows(
      rt_data %>%
        group_by(on_day) %>%
        summarize(n = mean(n)) %>%
        ungroup() %>%
        mutate(simulation = 'average')
    ) %>%
    mutate(type = ifelse(simulation == "average", "average", "simulations")) %>%
    ungroup() %>%
    filter(on_day <= 0.95*max(on_day)) %>%
    ggplot(aes(x = on_day, y = n, group = simulation, color = type, alpha = type)) +
    geom_line() + 
    scale_alpha_manual(values = c(1, 0.5)) +
    scale_color_manual(values = c('black', 'gray80')) +
    theme_bw() +
    xlab('day') +
    ggtitle("R_t by day") +
    geom_hline(yintercept = 1, linetype = 'dashed')
  
  list(
    sir_data = sir_data,
    rt_data = rt_data,
    sir_plot = sir_plot,
    i_plot = i_plot,
    rt_plot = rt_plot
  )
}

plot_single_simulation_vs_baseline <- function(
  simulation_result, benchmark_result
) {
  #' compare two simulation runs
}

plot_multiple_simulations_vs_baselines <- function(
  simulation_results, 
  benchmark_results
) {
  #' compare multiple simulation runs vs benchmark
}