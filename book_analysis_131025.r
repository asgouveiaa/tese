require(tidyverse)
require(deSolve)
require(ggpubr)

setwd("~/Documents/git/spillover25/.")
source("./new/prep_functions.R")

# Simulation parameters
n_years <- 10
total_days <- n_years * 365
tempo <- seq(1, total_days, by = 1)
lanina_years <- c(0)

precipitation_data <- generate_precipitation_dataset(
  years = n_years,
  la_nina_years = lanina_years,  
  nina_intensity = 1.5,
  drought_periods = create_periods(c(0,0,0), c(0,0,0)),
  flood_periods = create_periods(c(0,0), c(0,0)),
  flood_intensity = 3
) %>%
  mutate(p_normalizado = (precipitation - min(precipitation)) / 
           (max(precipitation) - min(precipitation)))

# PARÂMETROS PADRONIZADOS PARA TODAS AS COMBINAÇÕES
STANDARD_h <- 0.33
STANDARD_A_initial <- 40000
STANDARD_Ia_initial <- 4000
STANDARD_K_max <- 50000

# Configuração de cenários por combinação animal-patógeno
pathogen_config <- list(
  alouatta_mayaro = list(
    name = "Alouatta - Mayaro",
    animal = "A. seniculus",
    pathogen = "Mayaro",
    # Parâmetros específicos - MECANISMO: Transmissão vetorial
    encounters_day = 10, 
    bite_risk = 0.49, 
    p_V = 0.15, 
    mu = 0.0001, 
    r = 0.001,
    # OUTROS PARÂMETROS - USANDO VALORES PADRONIZADOS
    injury_risk = 0, w_injuries = 0, p_I = 0,
    eating_risk = 0, w_eating = 0, p_E = 0
  )
)

# Configuração de cenários epidemiológicos
scenario_config <- list(
  infection_hunting = list(
    name = "Com Infecção, Com Caça",
    Ia_initial_multiplier = 1,
    h_multiplier = 1
  )
)

# FUNÇÃO MODIFICADA PARA VARIAR DIAS DE CAÇA (INCLUINDO CAÇA A CADA 15 DIAS)
run_hunt_days_scenarios <- function() {
  
  # Definir cenários de dias de caça
  hunt_days_scenarios <- list(
    dois_dias = c(6, 7),           # Sábado e Domingo
    tre_dias = c(5, 6, 7), # Sexta, Sábado, Domingo
    cinco_dias = c(1:5),
    todos_dias = 1:7,               # Todos os dias da semana
    cada_15_dias = "every_15_days"  # Nova modalidade - caça a cada 15 dias
  )
  
  all_scenario_results <- list()
  counter <- 1
  
  # Executar para cada cenário de dias de caça
  for (hunt_days_name in names(hunt_days_scenarios)) {
    
    if (hunt_days_name == "cada_15_dias") {
      cat("Executando cenário: cada_15_dias - Caça a cada 15 dias\n")
      result <- run_every_15_days_scenario(
        pathogen_type = "alouatta_mayaro",
        scenario_type = "infection_hunting"
      )
    } else {
      cat("Executando cenário:", hunt_days_name, "- Dias:", 
          paste(hunt_days_scenarios[[hunt_days_name]], collapse = ","), "\n")
      
      result <- run_single_hunt_days_scenario(
        pathogen_type = "alouatta_mayaro",
        scenario_type = "infection_hunting",
        hunt_days_custom = hunt_days_scenarios[[hunt_days_name]]
      )
    }
    
    result$hunt_days_scenario <- hunt_days_name
    result$hunt_days_used <- ifelse(hunt_days_name == "cada_15_dias", 
                                    "every_15_days", 
                                    paste(hunt_days_scenarios[[hunt_days_name]], collapse = ","))
    result$num_hunt_days <- ifelse(hunt_days_name == "cada_15_dias", 
                                   round(365/15), # Aproximadamente 24 dias por ano
                                   length(hunt_days_scenarios[[hunt_days_name]]))
    
    all_scenario_results[[counter]] <- result
    counter <- counter + 1
  }
  
  # Combinar todos os resultados
  final_results <- bind_rows(all_scenario_results)
  
  cat("Total de cenários executados:", length(hunt_days_scenarios), "\n")
  
  return(final_results)
}

# FUNÇÃO PARA CENÁRIO DE CAÇA A CADA 15 DIAS
run_every_15_days_scenario <- function(pathogen_type, scenario_type) {
  
  pathogen_conf <- pathogen_config[[pathogen_type]]
  scenario_conf <- scenario_config[[scenario_type]]
  
  # USANDO VALORES PADRONIZADOS
  Ia_initial <- STANDARD_Ia_initial * scenario_conf$Ia_initial_multiplier
  h_value <- STANDARD_h * scenario_conf$h_multiplier
  
  # Parâmetros base - para caça a cada 15 dias usamos uma lógica diferente
  pars <- list(
    # HUNTING - Para caça a cada 15 dias, usamos hunt_days vazio e controlamos por outro mecanismo
    h = h_value,
    h_sazo = 1.2,
    base_time_forest = 1,
    time_max_optimal = 1,
    hunt_days = numeric(0),  # Vazio - controlaremos via lógica customizada
    hunter_experience = 0,
    hunt_frequency = 15,  # Novo parâmetro: frequência em dias
    
    # ANIMALS - PADRONIZADO
    r = pathogen_conf$r,
    K_max = STANDARD_K_max,
    K_sazo = 0.90,
    mu = pathogen_conf$mu,
    lambda_base = 0.1,
    lambda_agg_boost = 2.5,
    migration_rate = 0,
    
    # PATHOGENS
    phi = 50,
    nu = 0.7,
    P_max = 1e6,
    pathogen_seasonality = 0.3,
    
    # SPILLOVER - Direct transmission
    injurie_risk = pathogen_conf$injury_risk,
    p_I = pathogen_conf$p_I,
    w_injuries = pathogen_conf$w_injuries,
    
    processing_risk = 0, 
    p_P = 0, 
    w_processing = 0,
    
    eating_risk = pathogen_conf$eating_risk,
    p_E = pathogen_conf$p_E,
    w_eating = pathogen_conf$w_eating,
    
    cont_rate = 0, 
    p_F = 0, 
    w_fomites = 0,
    protective_equipment = 0,
    
    ## for consumers
    sharing_fraction = 0.75,
    home_processing_risk = 0,
    home_eating_risk = 0.1,
    home_ppe_effect = 0,
    cooking_protection = 0,
    
    # SPILLOVER - Vector transmission
    encounters_day = pathogen_conf$encounters_day,
    p_V = pathogen_conf$p_V,
    bite_risk = pathogen_conf$bite_risk,
    vector_seasonality = 0
  )
  
  # Condições iniciais PADRONIZADAS
  state <- c(
    Sh = 100,
    Ih = 0,
    Sc = 400,
    Ic = 0,
    A = STANDARD_A_initial,
    Ia = Ia_initial,
    P = 0
  )
  
  # Executar simulação com função modificada para caça a cada 15 dias
  out <- ode(y = state, times = tempo, func = spillover_model_every_15_days, 
             parms = pars, p_data = precipitation_data$p_normalizado)
  
  # Processar resultados
  resultado <- as.data.frame(out) %>%
    left_join(precipitation_data, by = c("time" = "day")) %>%
    mutate(
      pathogen_type = pathogen_type,
      scenario_type = scenario_type,
      animal = pathogen_conf$animal,
      pathogen = pathogen_conf$pathogen,
      scenario_name = scenario_conf$name,
      prevalence = ifelse((A + Ia) > 0, Ia / (A + Ia), 0),
      total_animals = A + Ia,
      has_infection = scenario_conf$Ia_initial_multiplier > 0,
      has_hunting = scenario_conf$h_multiplier > 0,
      # Informações sobre dias de caça
      hunt_days_used = "every_15_days",
      num_hunt_days = round(365/15),  # Aproximadamente 24 dias por ano
      h_used = h_value,
      A_initial_used = STANDARD_A_initial,
      Ia_initial_used = Ia_initial,
      K_max_used = STANDARD_K_max
    )
  
  return(resultado)
}

# FUNÇÃO DO MODELO MODIFICADA PARA CAÇA A CADA 15 DIAS
spillover_model_every_15_days <- function(t, state, pars, p_data) {
  with(as.list(c(state, pars)), {
    
    # Time indexing
    t_int <- pmax(1, pmin(floor(t), length(p_data)))
    p_norm <- p_data[t_int]
    
    # Climate-dependent effects with more complexity
    threshold <- 0.7
    clim_effect <- pmax(0, 0.3 * (p_norm - threshold) / (1 - threshold))
    agg_factor <- 1 + clim_effect
    
    # carrying capacity
    K_effective <- K_max * (K_sazo + (1 - K_sazo) * (1 - p_norm))
    
    lambda_effective <- lambda_base * (1 + (agg_factor - 1) * lambda_agg_boost)
    
    # LÓGICA DE CAÇA A CADA 15 DIAS
    # Determinar se é dia de caça (a cada 15 dias)
    is_hunting_day <- ifelse(floor(t) %% hunt_frequency == 0, 1, 0)
    
    # Weather-dependent hunting (hunters avoid rain)
    weather_factor <- ifelse(p_norm > 0.8, 0, 1.0)
    
    # More complex time in forest calculation
    optimal_climate <- 0.2
    climate_hunting_factor <- 1 - 0.5 * abs(p_norm - optimal_climate)
    
    time_in_forest_base <- (base_time_forest + 
                              (time_max_optimal - base_time_forest) * 
                              climate_hunting_factor) / 24
    
    # Effective rates on hunting days only
    time_in_forest <- is_hunting_day * time_in_forest_base * weather_factor
    hunting_rate <- is_hunting_day * h * (1 + h_sazo * p_norm)
    
    # Inicializa derivadas
    dSh <- dIh <- dSc <- dIc <- dA <- dIa <- dP <- 0
    hunter_spillover_prob <- 0
    consumer_spillover_prob <- 0
    
    # ------------------------------------------------------------
    # DINÂMICA DOS CAÇADORES (PRODUTORES)
    # ------------------------------------------------------------
    if (is_hunting_day == 1 && time_in_forest > 0) {
      animal_prop <- ifelse((Ia + A) > 0, Ia / (Ia + A), 0)
      
      # Componentes de risco para caçadores
      injury_risk_component <- injurie_risk * w_injuries * p_I * agg_factor * (1 - protective_equipment)
      processing_risk_hunter <- processing_risk * w_processing * p_P * (1 - protective_equipment)
      eating_risk_hunter <- eating_risk * w_eating * p_E
      fomite_risk <- time_in_forest * cont_rate * p_F * w_fomites * P * (1 - protective_equipment)
      
      # Risco total caçadores
      hunter_risk <- time_in_forest * hunting_rate * (
        injury_risk_component + 
          processing_risk_hunter + 
          eating_risk_hunter
      ) * animal_prop
      
      # Transmissão vetorial
      vector_seasonal_factor <- 1 + vector_seasonality * sin(2 * pi * (t - 120) / 365)
      vector_risk_hunter <- time_in_forest * encounters_day * 
        bite_risk * p_V * animal_prop * vector_seasonal_factor
      
      hunter_spillover <- hunter_risk * Sh
      vector_spillover <- vector_risk_hunter * Sh
      
      dSh <-  -(hunter_spillover + vector_spillover)
      dIh <-  hunter_spillover + vector_spillover
    }
    
    # ------------------------------------------------------------
    # DINÂMICA DOS FAMILIARES (CONSUMIDORES)
    # ------------------------------------------------------------
    if (hunting_rate > 0) {
      # Proporção de carne infectada compartilhada
      prop_infected <- ifelse((Ia + A) > 0, Ia / (Ia + A), 0)
      shared_meat <- sharing_fraction * hunting_rate * prop_infected
      
      # Riscos específicos para consumidores
      home_processing_risk_component <- shared_meat * home_processing_risk * (1 - home_ppe_effect)
      home_eating_risk_component <- shared_meat * home_eating_risk * (1 - cooking_protection)
      
      consumer_risk <- home_processing_risk_component + home_eating_risk_component
      
      consumer_spillover <- consumer_risk * Sc
      
      dSc <-  - consumer_spillover
      dIc <-  consumer_spillover
    }
    
    # ------------------------------------------------------------
    # DINÂMICA ANIMAL E AMBIENTAL
    # ------------------------------------------------------------
    # Transmissão animal-animal
    animal_transmission <- ifelse((Ia + A) > 0, 
                                  lambda_effective * A * (Ia / (Ia + A)), 0)
    
    # Equações diferenciais
    dA <- r * A * (1 - (A + Ia) / K_effective) - Sh * hunting_rate * A/(Ia+A) - mu * A - 
      animal_transmission + 0.1 * Ia
    
    dIa <- animal_transmission - Sh * hunting_rate * Ia/(Ia+A) - (mu*5) * Ia - 0.1 * Ia
    
    dP <- phi * Ia * (1 - (P / P_max)) - nu * P
    
    # ------------------------------------------------------------
    # OUTPUTS
    # ------------------------------------------------------------
    list(
      c(dSh = dSh, dIh = dIh, dSc = dSc, dIc = dIc, 
        dA = dA, dIa = dIa, dP = dP),
      
      hunter_spillover = ifelse(exists("hunter_spillover"), hunter_spillover, 0),
      vector_spillover = ifelse(exists("vector_spillover"), vector_spillover, 0),
      consumer_spillover = ifelse(exists("consumer_spillover"), consumer_spillover, 0),
      prop_meat_infected = ifelse((Ia + A) > 0, Ia / (Ia + A), 0),
      hunting_rate = hunting_rate,
      K_effective = K_effective,
      lambda_effective = lambda_effective,
      time_in_forest = time_in_forest,
      is_hunting_day = is_hunting_day,
      weekday = (floor(t) %% 7) + 1,  # Manter para compatibilidade
      hunter_spillover_prob = hunter_spillover_prob,
      consumer_spillover_prob = consumer_spillover_prob,
      hunt_frequency_used = hunt_frequency  # Para identificar este cenário
    )
  })
}

# FUNÇÃO PARA EXECUTAR CENÁRIO COM DIAS DE CAÇA PERSONALIZADOS (MANTIDA ORIGINAL)
run_single_hunt_days_scenario <- function(pathogen_type, scenario_type, 
                                          hunt_days_custom = c(6,7)) {
  
  pathogen_conf <- pathogen_config[[pathogen_type]]
  scenario_conf <- scenario_config[[scenario_type]]
  
  # USANDO VALORES PADRONIZADOS
  Ia_initial <- STANDARD_Ia_initial * scenario_conf$Ia_initial_multiplier
  h_value <- STANDARD_h * scenario_conf$h_multiplier
  
  # Parâmetros base COM DIAS DE CAÇA PERSONALIZADOS
  pars <- list(
    # HUNTING - COM DIAS PERSONALIZADOS
    h = h_value,
    h_sazo = 1.2,
    base_time_forest = 1,
    time_max_optimal = 1,
    hunt_days = hunt_days_custom,  # DIAS PERSONALIZADOS AQUI
    hunter_experience = 0,
    
    # ANIMALS - PADRONIZADO
    r = pathogen_conf$r,
    K_max = STANDARD_K_max,
    K_sazo = 0.90,
    mu = pathogen_conf$mu,
    lambda_base = 0.1,
    lambda_agg_boost = 2.5,
    migration_rate = 0,
    
    # PATHOGENS
    phi = 50,
    nu = 0.7,
    P_max = 1e6,
    pathogen_seasonality = 0.3,
    
    # SPILLOVER - Direct transmission
    injurie_risk = pathogen_conf$injury_risk,
    p_I = pathogen_conf$p_I,
    w_injuries = pathogen_conf$w_injuries,
    
    processing_risk = 0, 
    p_P = 0, 
    w_processing = 0,
    
    eating_risk = pathogen_conf$eating_risk,
    p_E = pathogen_conf$p_E,
    w_eating = pathogen_conf$w_eating,
    
    cont_rate = 0, 
    p_F = 0, 
    w_fomites = 0,
    protective_equipment = 0,
    
    ## for consumers
    sharing_fraction = 0.75,
    home_processing_risk = 0,
    home_eating_risk = 0.1,
    home_ppe_effect = 0,
    cooking_protection = 0,
    
    # SPILLOVER - Vector transmission
    encounters_day = pathogen_conf$encounters_day,
    p_V = pathogen_conf$p_V,
    bite_risk = pathogen_conf$bite_risk,
    vector_seasonality = 0
  )
  
  # Condições iniciais PADRONIZADAS
  state <- c(
    Sh = 100,
    Ih = 0,
    Sc = 400,
    Ic = 0,
    A = STANDARD_A_initial,
    Ia = Ia_initial,
    P = 0
  )
  
  # Executar simulação
  out <- ode(y = state, times = tempo, func = spillover_model, 
             parms = pars, p_data = precipitation_data$p_normalizado)
  
  # Processar resultados
  resultado <- as.data.frame(out) %>%
    left_join(precipitation_data, by = c("time" = "day")) %>%
    mutate(
      pathogen_type = pathogen_type,
      scenario_type = scenario_type,
      animal = pathogen_conf$animal,
      pathogen = pathogen_conf$pathogen,
      scenario_name = scenario_conf$name,
      prevalence = ifelse((A + Ia) > 0, Ia / (A + Ia), 0),
      total_animals = A + Ia,
      has_infection = scenario_conf$Ia_initial_multiplier > 0,
      has_hunting = scenario_conf$h_multiplier > 0,
      # Informações sobre dias de caça
      hunt_days_used = paste(hunt_days_custom, collapse = ","),
      num_hunt_days = length(hunt_days_custom),
      h_used = h_value,
      A_initial_used = STANDARD_A_initial,
      Ia_initial_used = Ia_initial,
      K_max_used = STANDARD_K_max
    )
  
  return(resultado)
}

# EXECUTAR CENÁRIOS DE DIAS DE CAÇA (INCLUINDO CADA 15 DIAS)
cat("Iniciando simulação de cenários de dias de caça...\n")
hunt_days_results <- run_hunt_days_scenarios()

# SALVAR RESULTADOS
#write_csv(hunt_days_results, "resultados_dias_caca_mayaro_com_15_dias.csv")

cat("Simulação concluída! Resultados salvos em 'resultados_dias_caca_mayaro_com_15_dias.csv'\n")
hunt_days_results <- read.csv('./resultados_dias_caca_mayaro_com_15_dias.csv')

# ANÁLISE RÁPIDA COMPARATIVA
summary_hunt_days <- hunt_days_results %>%
  group_by(hunt_days_scenario, num_hunt_days, hunt_days_used) %>%
  summarise(
    max_cacadores_infectados = max(Ih),
    total_cacadores_infectados = max(cumsum(Ih)),
    prevalencia_maxima = max(prevalence, na.rm = TRUE),
    prevalencia_media = mean(prevalence, na.rm = TRUE),
    total_animais_final = last(total_animals),
    dias_ate_primeira_infeccao = ifelse(any(Ih > 0), time[which(Ih > 0)[1]], NA),
    total_dias_caca_efetivos = sum(hunting_rate > 0),
    .groups = 'drop'
  ) %>%
  arrange(num_hunt_days)

print(summary_hunt_days)

hunt_days_results %>%
  ggplot() +
  geom_line(aes(y=Ia/(Ia+A)*100, x=Ih/(Ih+Sh)*100, color = hunt_days_scenario))+
  labs(title = "",
       y = "Prevalence (%)", x = "Time (years)",
       color = "") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top") +
  facet_wrap(~hunt_days_scenario)
  # scale_x_continuous(
  #   breaks = seq(1, max(all_results$time), by = 365),  # Breaks at 0, 365, 730
  #   labels = c(seq(1,10,1))
  # )


vary_freq <- summary_hunt_days %>%
  mutate(hunt_days_scenario = factor(hunt_days_scenario,
                levels = c("cada_15_dias","dois_dias","tre_dias",
                           "cinco_dias","todos_dias"),
                labels = c("2x month", "2x week", "3x week", "5x week", "7x week"))) %>% 
  ggplot(aes(x=factor(hunt_days_scenario), y=max_cacadores_infectados/100, group = 1), color = "black") +
  geom_point() +
  geom_line() +
  theme_bw(base_size = 12) + 
  labs(x = "Hunting frequency", y = "Proportion of infected hunters") 


# varying number of hunters

# FUNÇÃO PRINCIPAL PARA HEATMAP - VARIA DIAS DE CAÇA E NÚMERO DE CAÇADORES
run_heatmap_scenarios <- function() {
  
  # Definir ranges para o heatmap
  sh_range <- seq(10, 500, by = 5)  # 10, 20, 30, ..., 100 caçadores
  hunt_days_scenarios <- list(
    dois_dias = c(6, 7),           # Sábado e Domingo
    tre_dias = c(5, 6, 7), # Sexta, Sábado, Domingo
    cinco_dias = c(1:5),
    todos_dias = 1:7,               # Todos os dias da semana
    cada_15_dias = "every_15_days"  # Nova modalidade - caça a cada 15 dias
  )
  
  all_scenario_results <- list()
  counter <- 1
  total_combinations <- length(sh_range) * length(hunt_days_scenarios)
  current_combination <- 0
  
  cat("Iniciando simulação para heatmap:", total_combinations, "combinações\n")
  
  # Loop através de todas as combinações
  for (sh_value in sh_range) {
    for (hunt_days_name in names(hunt_days_scenarios)) {
      
      current_combination <- current_combination + 1
      cat("Progresso:", current_combination, "/", total_combinations, 
          "- Sh:", sh_value, "- Dias:", hunt_days_name, "\n")
      
      if (hunt_days_name == "cada_15_dias") {
        result <- run_every_15_days_scenario(
          pathogen_type = "alouatta_mayaro",
          scenario_type = "infection_hunting",
          sh_initial = sh_value
        )
      } else {
        result <- run_single_hunt_days_scenario(
          pathogen_type = "alouatta_mayaro",
          scenario_type = "infection_hunting",
          hunt_days_custom = hunt_days_scenarios[[hunt_days_name]],
          sh_initial = sh_value
        )
      }
      
      # Adicionar metadados para o heatmap
      result$sh_initial <- sh_value
      result$hunt_days_scenario <- hunt_days_name
      result$hunt_days_used <- ifelse(hunt_days_name == "cada_15_dias", 
                                      "every_15_days", 
                                      paste(hunt_days_scenarios[[hunt_days_name]], collapse = ","))
      result$num_hunt_days <- ifelse(hunt_days_name == "cada_15_dias", 
                                     round(365/15),
                                     length(hunt_days_scenarios[[hunt_days_name]]))
      
      all_scenario_results[[counter]] <- result
      counter <- counter + 1
    }
  }
  
  # Combinar todos os resultados
  final_results <- bind_rows(all_scenario_results)
  
  cat("Simulação do heatmap concluída! Total de cenários:", nrow(final_results), "\n")
  
  return(final_results)
}

# FUNÇÃO MODIFICADA PARA CAÇA A CADA 15 DIAS (COM Sh VARIÁVEL)
run_every_15_days_scenario <- function(pathogen_type, scenario_type, sh_initial = 100) {
  
  pathogen_conf <- pathogen_config[[pathogen_type]]
  scenario_conf <- scenario_config[[scenario_type]]
  
  Ia_initial <- STANDARD_Ia_initial * scenario_conf$Ia_initial_multiplier
  h_value <- STANDARD_h * scenario_conf$h_multiplier
  
  pars <- list(
    h = h_value,
    h_sazo = 1.2,
    base_time_forest = 1,
    time_max_optimal = 1,
    hunt_days = numeric(0),
    hunter_experience = 0,
    hunt_frequency = 15,
    
    r = pathogen_conf$r,
    K_max = STANDARD_K_max,
    K_sazo = 0.90,
    mu = pathogen_conf$mu,
    lambda_base = 0.1,
    lambda_agg_boost = 2.5,
    migration_rate = 0,
    
    phi = 50,
    nu = 0.7,
    P_max = 1e6,
    pathogen_seasonality = 0.3,
    
    injurie_risk = pathogen_conf$injury_risk,
    p_I = pathogen_conf$p_I,
    w_injuries = pathogen_conf$w_injuries,
    
    processing_risk = 0, p_P = 0, w_processing = 0,
    eating_risk = pathogen_conf$eating_risk,
    p_E = pathogen_conf$p_E,
    w_eating = pathogen_conf$w_eating,
    
    cont_rate = 0, p_F = 0, w_fomites = 0,
    protective_equipment = 0,
    
    sharing_fraction = 0.75,
    home_processing_risk = 0,
    home_eating_risk = 0.1,
    home_ppe_effect = 0,
    cooking_protection = 0,
    
    encounters_day = pathogen_conf$encounters_day,
    p_V = pathogen_conf$p_V,
    bite_risk = pathogen_conf$bite_risk,
    vector_seasonality = 0
  )
  
  # Condições iniciais COM Sh VARIÁVEL
  state <- c(
    Sh = sh_initial,  # AGORA VARIÁVEL
    Ih = 0,
    Sc = 400,
    Ic = 0,
    A = STANDARD_A_initial,
    Ia = Ia_initial,
    P = 0
  )
  
  out <- ode(y = state, times = tempo, func = spillover_model_every_15_days, 
             parms = pars, p_data = precipitation_data$p_normalizado)
  
  resultado <- as.data.frame(out) %>%
    left_join(precipitation_data, by = c("time" = "day")) %>%
    mutate(
      pathogen_type = pathogen_type,
      scenario_type = scenario_type,
      animal = pathogen_conf$animal,
      pathogen = pathogen_conf$pathogen,
      scenario_name = scenario_conf$name,
      prevalence = ifelse((A + Ia) > 0, Ia / (A + Ia), 0),
      total_animals = A + Ia,
      has_infection = scenario_conf$Ia_initial_multiplier > 0,
      has_hunting = scenario_conf$h_multiplier > 0,
      sh_initial = sh_initial,  # Registrar o valor usado
      hunt_days_used = "every_15_days",
      num_hunt_days = round(365/15),
      h_used = h_value,
      A_initial_used = STANDARD_A_initial,
      Ia_initial_used = Ia_initial,
      K_max_used = STANDARD_K_max
    )
  
  return(resultado)
}

# FUNÇÃO MODIFICADA PARA DIAS DE CAÇA REGULARES (COM Sh VARIÁVEL)
run_single_hunt_days_scenario <- function(pathogen_type, scenario_type, 
                                          hunt_days_custom = c(6,7),
                                          sh_initial = 100) {  # AGORA COM Sh VARIÁVEL
  
  pathogen_conf <- pathogen_config[[pathogen_type]]
  scenario_conf <- scenario_config[[scenario_type]]
  
  Ia_initial <- STANDARD_Ia_initial * scenario_conf$Ia_initial_multiplier
  h_value <- STANDARD_h * scenario_conf$h_multiplier
  
  pars <- list(
    h = h_value,
    h_sazo = 1.2,
    base_time_forest = 1,
    time_max_optimal = 1,
    hunt_days = hunt_days_custom,
    hunter_experience = 0,
    
    r = pathogen_conf$r,
    K_max = STANDARD_K_max,
    K_sazo = 0.90,
    mu = pathogen_conf$mu,
    lambda_base = 0.1,
    lambda_agg_boost = 2.5,
    migration_rate = 0,
    
    phi = 50,
    nu = 0.7,
    P_max = 1e6,
    pathogen_seasonality = 0.3,
    
    injurie_risk = pathogen_conf$injury_risk,
    p_I = pathogen_conf$p_I,
    w_injuries = pathogen_conf$w_injuries,
    
    processing_risk = 0, p_P = 0, w_processing = 0,
    eating_risk = pathogen_conf$eating_risk,
    p_E = pathogen_conf$p_E,
    w_eating = pathogen_conf$w_eating,
    
    cont_rate = 0, p_F = 0, w_fomites = 0,
    protective_equipment = 0,
    
    sharing_fraction = 0.75,
    home_processing_risk = 0,
    home_eating_risk = 0.1,
    home_ppe_effect = 0,
    cooking_protection = 0,
    
    encounters_day = pathogen_conf$encounters_day,
    p_V = pathogen_conf$p_V,
    bite_risk = pathogen_conf$bite_risk,
    vector_seasonality = 0
  )
  
  # Condições iniciais COM Sh VARIÁVEL
  state <- c(
    Sh = sh_initial,  # AGORA VARIÁVEL
    Ih = 0,
    Sc = 400,
    Ic = 0,
    A = STANDARD_A_initial,
    Ia = Ia_initial,
    P = 0
  )
  
  out <- ode(y = state, times = tempo, func = spillover_model, 
             parms = pars, p_data = precipitation_data$p_normalizado)
  
  resultado <- as.data.frame(out) %>%
    left_join(precipitation_data, by = c("time" = "day")) %>%
    mutate(
      pathogen_type = pathogen_type,
      scenario_type = scenario_type,
      animal = pathogen_conf$animal,
      pathogen = pathogen_conf$pathogen,
      scenario_name = scenario_conf$name,
      prevalence = ifelse((A + Ia) > 0, Ia / (A + Ia), 0),
      total_animals = A + Ia,
      has_infection = scenario_conf$Ia_initial_multiplier > 0,
      has_hunting = scenario_conf$h_multiplier > 0,
      sh_initial = sh_initial,  # Registrar o valor usado
      hunt_days_used = paste(hunt_days_custom, collapse = ","),
      num_hunt_days = length(hunt_days_custom),
      h_used = h_value,
      A_initial_used = STANDARD_A_initial,
      Ia_initial_used = Ia_initial,
      K_max_used = STANDARD_K_max
    )
  
  return(resultado)
}

# EXECUTAR SIMULAÇÃO PARA HEATMAP
cat("Iniciando simulação para heatmap...\n")
heatmap_results <- run_heatmap_scenarios()

# SALVAR RESULTADOS COMPLETOS
#write_csv(heatmap_results, "dados_heatmap_completo.csv")

heatmap_results <- read.csv("./dados_heatmap_completo.csv")

heatmap_results2 <- heatmap_results %>%
  mutate(A = ifelse(A<0,0,A),
         Ia = ifelse(Ia<0,0,Ia),
         P = ifelse(P<0,0,P))

vary_freq_sh <- heatmap_results2 %>% 
  mutate(hunt_days_scenario = factor(hunt_days_scenario,
                                     levels = c("cada_15_dias","dois_dias","tre_dias",
                                                "cinco_dias","todos_dias"),
                                     labels = c("2x month", "2x week", "3x week", "5x week", "7x week"))) %>% 
  group_by(sh_initial, hunt_days_scenario) %>% 
  summarise(final_Ih = tail(Ih,1),
            final_Sh = tail(Sh,1)) %>%
  mutate(prev = final_Ih/(final_Ih + final_Sh)) %>% 
  ggplot() + 
  geom_tile(aes(x=hunt_days_scenario, y = factor(sh_initial), fill = prev)) + 
  scale_fill_viridis_c(limits = c(0, 1), 
                       option = "mako", 
                       breaks = c(0.25, 0.50, 0.75),
                       labels = c("0.25", "0.50", "0.75"))+
  theme_bw(base_size = 12) + 
  theme(legend.position = "top") +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  scale_y_discrete(
    expand = c(0,0),
    breaks = levels(factor(heatmap_results2$sh_initial))[seq(1, 
                                                             length(levels(factor(heatmap_results2$sh_initial))), 
                                                             by = 15)]  # A cada 2 valores
  ) + 
  labs(x = "Hunter Frequency",
       y = "Hunter population",
       fill = "Proportion of infected hunters" ); vary_freq_sh# +
  # theme(legend.position = "right")


scenario1 <- ggarrange(vary_freq,vary_freq_sh, align = "hv", labels = "auto")
ggsave("./book_chapter/figs/fig.scen1.png",dpi = 600)

#variando h e lambda

# FUNÇÃO PARA VARIAR lambda_base
run_lambda_base_scenarios <- function() {
  
  # Valores de lambda_base para testar
  lambda_base_range <- seq(0, 1, length.out=100)

  
  all_scenario_results <- list()
  counter <- 1
  
  cat("Cenário 1: Variando lambda_base\n")
  
  for (lambda_val in lambda_base_range) {
    cat("  lambda_base =", lambda_val, "\n")
    
    result <- run_single_parameter_scenario(
      pathogen_type = "alouatta_mayaro",
      scenario_type = "infection_hunting",
      lambda_base_custom = lambda_val,
      h_custom = 0.33  # Valor padrão
    )
    
    result$scenario_group <- "lambda_base"
    result$parameter_name <- "lambda_base"
    result$parameter_value <- lambda_val
    result$h_value <- 0.33
    
    all_scenario_results[[counter]] <- result
    counter <- counter + 1
  }
  
  return(bind_rows(all_scenario_results))
}

# FUNÇÃO PARA VARIAR lambda_base E h CONJUNTAMENTE
run_lambda_h_scenarios <- function() {
  
  # Valores para testar
  lambda_base_range <- seq(0,0.5,length.out=20)
  h_range <- seq(0,0.6,length.out=20)
  
  all_scenario_results <- list()
  counter <- 1
  total_combinations <- length(lambda_base_range) * length(h_range)
  current_combination <- 0
  
  cat("Cenário 2: Variando lambda_base × h\n")
  cat("Total de combinações:", total_combinations, "\n")
  
  for (lambda_val in lambda_base_range) {
    for (h_val in h_range) {
      
      current_combination <- current_combination + 1
      cat("  Progresso:", current_combination, "/", total_combinations,
          "- lambda_base:", lambda_val, "- h:", h_val, "\n")
      
      result <- run_single_parameter_scenario(
        pathogen_type = "alouatta_mayaro",
        scenario_type = "infection_hunting",
        lambda_base_custom = lambda_val,
        h_custom = h_val
      )
      
      result$scenario_group <- "lambda_h_combination"
      result$parameter_name <- paste0("lambda_", lambda_val, "_h_", h_val)
      result$lambda_base_value <- lambda_val
      result$h_value <- h_val
      
      all_scenario_results[[counter]] <- result
      counter <- counter + 1
    }
  }
  
  return(bind_rows(all_scenario_results))
}

# FUNÇÃO PARA EXECUTAR CENÁRIO COM PARÂMETROS PERSONALIZADOS
run_single_parameter_scenario <- function(pathogen_type, scenario_type,
                                          lambda_base_custom = 0.1,
                                          h_custom = 0.33) {
  
  pathogen_conf <- pathogen_config[[pathogen_type]]
  scenario_conf <- scenario_config[[scenario_type]]
  
  Ia_initial <- STANDARD_Ia_initial * scenario_conf$Ia_initial_multiplier
  h_value <- h_custom * scenario_conf$h_multiplier  # Usar valor customizado
  
  # Parâmetros base COM PARÂMETROS PERSONALIZADOS
  pars <- list(
    # HUNTING - COM h PERSONALIZADO
    h = h_value,
    h_sazo = 1.2,
    base_time_forest = 1,
    time_max_optimal = 1,
    hunt_days = c(6, 7),  # Padrão: fim de semana
    hunter_experience = 0,
    
    # ANIMALS - COM lambda_base PERSONALIZADO
    r = pathogen_conf$r,
    K_max = STANDARD_K_max,
    K_sazo = 0.90,
    mu = pathogen_conf$mu,
    lambda_base = lambda_base_custom,  # PARÂMETRO PERSONALIZADO
    lambda_agg_boost = 2.5,
    migration_rate = 0,
    
    # PATHOGENS
    phi = 50,
    nu = 0.7,
    P_max = 1e6,
    pathogen_seasonality = 0.3,
    
    # SPILLOVER - Direct transmission
    injurie_risk = pathogen_conf$injury_risk,
    p_I = pathogen_conf$p_I,
    w_injuries = pathogen_conf$w_injuries,
    
    processing_risk = 0, 
    p_P = 0, 
    w_processing = 0,
    
    eating_risk = pathogen_conf$eating_risk,
    p_E = pathogen_conf$p_E,
    w_eating = pathogen_conf$w_eating,
    
    cont_rate = 0, 
    p_F = 0, 
    w_fomites = 0,
    protective_equipment = 0,
    
    ## for consumers
    sharing_fraction = 0.75,
    home_processing_risk = 0,
    home_eating_risk = 0.1,
    home_ppe_effect = 0,
    cooking_protection = 0,
    
    # SPILLOVER - Vector transmission
    encounters_day = pathogen_conf$encounters_day,
    p_V = pathogen_conf$p_V,
    bite_risk = pathogen_conf$bite_risk,
    vector_seasonality = 0
  )
  
  # Condições iniciais PADRONIZADAS
  state <- c(
    Sh = 100,
    Ih = 0,
    Sc = 400,
    Ic = 0,
    A = STANDARD_A_initial,
    Ia = Ia_initial,
    P = 0
  )
  
  # Executar simulação
  out <- ode(y = state, times = tempo, func = spillover_model, 
             parms = pars, p_data = precipitation_data$p_normalizado)
  
  # Processar resultados
  resultado <- as.data.frame(out) %>%
    left_join(precipitation_data, by = c("time" = "day")) %>%
    mutate(
      pathogen_type = pathogen_type,
      scenario_type = scenario_type,
      animal = pathogen_conf$animal,
      pathogen = pathogen_conf$pathogen,
      scenario_name = scenario_conf$name,
      prevalence = ifelse((A + Ia) > 0, Ia / (A + Ia), 0),
      total_animals = A + Ia,
      has_infection = scenario_conf$Ia_initial_multiplier > 0,
      has_hunting = scenario_conf$h_multiplier > 0,
      # Informações sobre parâmetros usados
      lambda_base_used = lambda_base_custom,
      h_used = h_value,
      A_initial_used = STANDARD_A_initial,
      Ia_initial_used = Ia_initial,
      K_max_used = STANDARD_K_max
    )
  
  return(resultado)
}

# FUNÇÃO DO MODELO (ORIGINAL)
spillover_model <- function(t, state, pars, p_data) {
  with(as.list(c(state, pars)), {
    
    # Time indexing
    t_int <- pmax(1, pmin(floor(t), length(p_data)))
    p_norm <- p_data[t_int]
    
    # Climate-dependent effects with more complexity
    threshold <- 0.7
    clim_effect <- pmax(0, 0.3 * (p_norm - threshold) / (1 - threshold))
    agg_factor <- 1 + clim_effect
    
    # carrying capacity
    K_effective <- K_max * (K_sazo + (1 - K_sazo) * (1 - p_norm))
    
    lambda_effective <- lambda_base * (1 + (agg_factor - 1) * lambda_agg_boost)
    
    # Hunting schedule with experience factor
    weekday <- (floor(t) %% 7) + 1  
    is_hunting_day <- ifelse(weekday %in% hunt_days, 1, 0)
    
    # Weather-dependent hunting (hunters avoid rain)
    weather_factor <- ifelse(p_norm > 0.8, 0, 1.0)
    
    # More complex time in forest calculation
    optimal_climate <- 0.2
    climate_hunting_factor <- 1 - 0.5 * abs(p_norm - optimal_climate)
    
    time_in_forest_base <- (base_time_forest + 
                              (time_max_optimal - base_time_forest) * 
                              climate_hunting_factor) / 24
    
    # Effective rates on hunting days only
    time_in_forest <- is_hunting_day * time_in_forest_base * weather_factor
    hunting_rate <- is_hunting_day * h * (1 + h_sazo * p_norm)
    
    # Inicializa derivadas
    dSh <- dIh <- dSc <- dIc <- dA <- dIa <- dP <- 0
    hunter_spillover_prob <- 0
    consumer_spillover_prob <- 0
    
    # ------------------------------------------------------------
    # DINÂMICA DOS CAÇADORES (PRODUTORES)
    # ------------------------------------------------------------
    if (is_hunting_day == 1 && time_in_forest > 0) {
      animal_prop <- ifelse((Ia + A) > 0, Ia / (Ia + A), 0)
      
      # Componentes de risco para caçadores
      injury_risk_component <- injurie_risk * w_injuries * p_I * agg_factor * (1 - protective_equipment)
      processing_risk_hunter <- processing_risk * w_processing * p_P * (1 - protective_equipment)
      eating_risk_hunter <- eating_risk * w_eating * p_E
      fomite_risk <- time_in_forest * cont_rate * p_F * w_fomites * P * (1 - protective_equipment)
      
      # Risco total caçadores
      hunter_risk <- time_in_forest * hunting_rate * (
        injury_risk_component + 
          processing_risk_hunter + 
          eating_risk_hunter
      ) * animal_prop
      
      # Transmissão vetorial
      vector_seasonal_factor <- 1 + vector_seasonality * sin(2 * pi * (t - 120) / 365)
      vector_risk_hunter <- time_in_forest * encounters_day * 
        bite_risk * p_V * animal_prop * vector_seasonal_factor
      
      hunter_spillover <- hunter_risk * Sh
      vector_spillover <- vector_risk_hunter * Sh
      
      dSh <-  -(hunter_spillover + vector_spillover)
      dIh <-  hunter_spillover + vector_spillover
    }
    
    # ------------------------------------------------------------
    # DINÂMICA DOS FAMILIARES (CONSUMIDORES)
    # ------------------------------------------------------------
    if (hunting_rate > 0) {
      # Proporção de carne infectada compartilhada
      prop_infected <- ifelse((Ia + A) > 0, Ia / (Ia + A), 0)
      shared_meat <- sharing_fraction * hunting_rate * prop_infected
      
      # Riscos específicos para consumidores
      home_processing_risk_component <- shared_meat * home_processing_risk * (1 - home_ppe_effect)
      home_eating_risk_component <- shared_meat * home_eating_risk * (1 - cooking_protection)
      
      consumer_risk <- home_processing_risk_component + home_eating_risk_component
      
      consumer_spillover <- consumer_risk * Sc
      
      dSc <-  - consumer_spillover
      dIc <-  consumer_spillover
    }
    
    # ------------------------------------------------------------
    # DINÂMICA ANIMAL E AMBIENTAL
    # ------------------------------------------------------------
    # Transmissão animal-animal
    animal_transmission <- ifelse((Ia + A) > 0, 
                                  lambda_effective * A * (Ia / (Ia + A)), 0)
    
    # Equações diferenciais
    dA <- r * A * (1 - (A + Ia) / K_effective) - Sh * hunting_rate * A/(Ia+A) - mu * A - 
      animal_transmission + 0.1 * Ia
    
    dIa <- animal_transmission - Sh * hunting_rate * Ia/(Ia+A) - (mu*5) * Ia - 0.1 * Ia
    
    dP <- phi * Ia * (1 - (P / P_max)) - nu * P
    
    # ------------------------------------------------------------
    # OUTPUTS
    # ------------------------------------------------------------
    list(
      c(dSh = dSh, dIh = dIh, dSc = dSc, dIc = dIc, 
        dA = dA, dIa = dIa, dP = dP),
      
      hunter_spillover = ifelse(exists("hunter_spillover"), hunter_spillover, 0),
      vector_spillover = ifelse(exists("vector_spillover"), vector_spillover, 0),
      consumer_spillover = ifelse(exists("consumer_spillover"), consumer_spillover, 0),
      prop_meat_infected = ifelse((Ia + A) > 0, Ia / (Ia + A), 0),
      hunting_rate = hunting_rate,
      K_effective = K_effective,
      lambda_effective = lambda_effective,
      time_in_forest = time_in_forest,
      is_hunting_day = is_hunting_day,
      weekday = weekday,
      hunter_spillover_prob = hunter_spillover_prob,
      consumer_spillover_prob = consumer_spillover_prob
    )
  })
}

# EXECUTAR OS DOIS CENÁRIOS
cat("Iniciando simulações para análise de parâmetros...\n")

# Cenário 1: Variar apenas lambda_base
cat("\n=== CENÁRIO 1: Variando lambda_base ===\n")
lambda_results <- run_lambda_base_scenarios()

# Cenário 2: Variar lambda_base e h conjuntamente
cat("\n=== CENÁRIO 2: Variando lambda_base × h ===\n")
lambda_h_results <- run_lambda_h_scenarios()

# Combinar resultados
all_parameter_results <- bind_rows(
  lambda_results,
  lambda_h_results
)

# SALVAR RESULTADOS
write_csv(all_parameter_results, "resultados_parametros_lambda_h.csv")

# CRIAR RESUMOS PARA ANÁLISE
summary_lambda <- lambda_results %>%
  group_by(lambda_base_used) %>%
  summarise(
    prevalencia_maxima = max(prevalence, na.rm = TRUE),
    prevalencia_media = mean(prevalence, na.rm = TRUE),
    total_cacadores_infectados = max(Ih),
    pico_prevalencia = time[which.max(prevalence)],
    total_animais_final = last(total_animals),
    .groups = 'drop'
  )

summary_lambda_h <- lambda_h_results %>%
  group_by(lambda_base_value, h_value) %>%
  summarise(
    prevalencia_maxima = max(prevalence, na.rm = TRUE),
    total_cacadores_infectados = max(Ih),
    taxa_ataque = max(Ih) / 100,  # Assume 100 caçadores iniciais
    estabilidade_prevalencia = sd(prevalence[time > 1000], na.rm = TRUE),
    .groups = 'drop'
  )

# SALVAR RESUMOS
write_csv(summary_lambda, "resumo_lambda_base.csv")
write_csv(summary_lambda_h, "resumo_lambda_h.csv")

# VISUALIZAÇÕES
library(viridis)


lambda_results %>% 
  ggplot(aes(y=Ia,x=Ih, color = lambda_base_used)) + 
  geom_line()

# Gráfico 1: Efeito do lambda_base na prevalência
ggplot(summary_lambda, aes(x = lambda_base_used, y = prevalencia_maxima)) +
  # geom_point(size = 3, color = "#2E86AB") +
  geom_point(col = "black") +
  geom_line(col = "black")+
  labs(title = "Efeito da Taxa de Transmissão (lambda_base) na Prevalência Máxima",
       x = "lambda_base (taxa de transmissão)",
       y = "Prevalência Máxima") +
  theme_minimal()



 # Heatmap: lambda_base × h
p2 <- ggplot(summary_lambda_h, aes(x = factor(lambda_base_value),
                                   y = factor(h_value),
                                   fill = prevalencia_maxima)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako"); p2
#   geom_text(aes(label = round(prevalencia_maxima, 3)))

lambda_h_results %>% 
  mutate(A = ifelse(A<0,0,A),
         Ia = ifelse(Ia<0,0,Ia),
         P = ifelse(P<0,0,P)) %>% 
  group_by(h_value, lambda_base_value) %>% 
  summarise(final_Ih = tail(Ih,1),
            final_Sh = tail(Sh,1)) %>%
  mutate(prev = final_Ih/(final_Ih + final_Sh)*100) %>% 
  ggplot() +
  geom_tile(aes(x=factor(lambda_base_value),y=factor(h_value), fill = prev)) +
  scale_fill_viridis_c(option = "mako")

fig.lambda_h <- lambda_h_results %>% 
  ggplot() +
  geom_tile(aes(x=lambda_base_value,y=h_value, fill = prevalence)) + 
  scale_fill_viridis_c(option = "mako")

ggsave("./fig.lambda_h.png",fig.lambda_h)



t <- read.csv2("~/Downloads/sensity2parms.csv")

scen2_p2 <- ggplot(t) +
  geom_tile(aes(y=param2_value, x =param1_value, fill = final_prevalenceH)) +
  scale_fill_viridis_c(, option = "mako") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "top") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "λ0", y = "h0", fill = "Proportion of infected hunters")

scen2_p1 <- ggplot(summary_lambda, aes(x = lambda_base_used, y = prevalencia_maxima)) +
  # geom_point(size = 3, color = "#2E86AB") +
  geom_point(col = "black") +
  geom_line(col = "black")+
  labs(title = "",
       x = "λ0",
       y = "Proportion of infected hunters") +
  theme_bw(base_size = 12); scen2_p1



scenario2 <- ggarrange(scen2_p1,scen2_p2, align = "hv", labels = "auto")
ggsave("./book_chapter/figs/fig.scen2.png",scenario2,dpi = 900)
