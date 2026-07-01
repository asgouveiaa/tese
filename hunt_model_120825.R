#====PACOTES E FUNCOES EXTRAS====
pacman::p_load(
    tidyverse,deSolve,
    lubridate,magrittr,
    ggpubr
)

source("./hunt/prep_functions.R")

#====PARAMETROS E CONDICOES INICIAIS====
n_years <- 2
total_days <- n_years * 365
tempo <- seq(1, total_days, by = 1)
lanina_years <- c(0)
optimal_climate = 0.4

pars <- list(  
  # HUNTING
  h = 0.05,                  # hunting efficiency (scenario-dependent)
  h_sazo = 0.1,             # seasonal increase 
  base_time_forest = 4,      # average daily hunting time (hours)
  time_max_optimal = 6,      # increased max hunting time
  hunt_days = c(6,7),        # weekend hunting 
  hunter_experience = 0,     # experience factor
 
  
  # ANIMALS
  r = 0.25,                  #10seasonal decrease,
  K_max = 1500,              # maximum carrying capacity
  K_sazo = 0.85,             # seasonal decrease
  mu = 0.05,                 # daily mortality
  lambda_base = 0.08,        # baseline transmission rate
  lambda_agg_boost = 1.2,    # stronger aggregation effect
  #migration_rate = 0,        # migration rate  
  
  # PATHOGENS
  phi = 50,                   # shedding rate
  nu = 0.7,                   # slower pathogen decay
  P_max = 1e6,                # max pathogen load
  pathogen_seasonality = 0.3, #seasonal pathogen survival
  
  #HUNTERS (DIRECT TRANSMISSION)
  injury_risk = 0, prob_trans_injury = 0, weight_injury = 0,        #  injury risk
  processing_risk = 0, prob_trans_processing = 0, weight_processing = 0,   #  processing risk
  eating_risk = 0.5, prob_trans_eating = 1, weight_eating = 1,         #  eating risk
  cont_rate = 0, prob_trans_fomite = 0, weight_fomites = 0,#,           #  fomite risk
  #protective_equipment = 0,                        #  PPE usage rate
 
  #HUNTER (VECTOR-BORNE)
  encounters_day = 0,      #  number of vector encounter per day
  prob_trans_vector = 0.05,               #  higher infection probability
  bite_risk = 0.2,         #  higher bite risk
  #vector_seasonality = 0   #  seasonal vector activity
  
  #CONSUMERS
  sharing_fraction = 0.75,      # meat shared fraction
  home_processing_risk = 0,     # home processing risk
  home_eating_risk = 0.25,      # home eating risk
  home_ppe_effect = 0,          # PPE at home
  cooking_protection = 0       # risk reduction by cooking
)

# Initial conditions
state <- c(
  Sh = 100,               # Hunters S
  Ih = 0,                 # Hunters I
  Sc = 400,               # Consumers S
  Ic = 0,                 # Consumers I
  Sa = 800,                # Animals S
  Ia = 100,               # Animals I
  P = 0                   # Pathogen in the enviroment
)

#====GERANDO PRECIPITACAO====
precipitation_data <- generate_precipitation_dataset(
  years = n_years,
  la_nina_years = lanina_years,  
  nina_intensity = 1.5,
  drought_periods = create_periods(c(0,0,0), c(0,0,0)),
  flood_periods = create_periods(c(0,0), c(0,0)),
  flood_intensity = 3
) |>
  mutate(p_normalizado =
   (precipitation - min(precipitation)) /
    (max(precipitation) - min(precipitation)))

#====DEFINICAO MODELO====
spillover_model <- function(t, state, pars, p_data) {
  with(as.list(c(state, pars)), {
    
    # Time indexing
    t_int <- pmax(1, pmin(floor(t), length(p_data)))
    p_norm <- p_data[t_int]
    
    # Climate-dependent effects with more complexity
    # Unified climate effect (0-1 scale where 0=unfavorable, 1=favorable)
    climate_effect <- 1 - abs(p_norm - optimal_climate)  # optimal_climate could be 0.4 as in your example
    
    # Agreggation factor
    agg_factor <- 1 + climate_effect


    # Apply this single effect to all climate-dependent components:
    # 1. Pathogen transmission
    lambda_effective <- lambda_base * (1 + (agg_factor - 1) * lambda_agg_boost) 

    # 2. Carrying capacity
    K_effective <- K_max * (K_sazo + (1 - K_sazo) * (1-climate_effect))

    # 3. Hunting activity
    climate_hunting_factor <- 1 - 0.5 * abs(p_norm - optimal_climate)
    weather_factor <- ifelse(p_norm > 0.8, 0, climate_effect)  # Still keep heavy rain cutoff
    time_in_forest_base <- (base_time_forest + 
                        (time_max_optimal - base_time_forest) * 
                        climate_hunting_factor) / 24
    
    # Hunting schedule with experience factor
    weekday <- (floor(t) %% 7) + 1  
    is_hunting_day <- ifelse(weekday %in% hunt_days, 1, 0)
    
    # Effective rates on hunting days only
    time_in_forest <- is_hunting_day * time_in_forest_base * weather_factor
    hunting_rate <- is_hunting_day * h * (1 + h_sazo * p_norm) #* hunter_experience
    
    # Start derivatives
    dSh <- dIh <- dSc <- dIc <- dSa <- dIa <- dP <- 0
    hunter_spillover_prob <- 0
    consumer_spillover_prob <- 0

      
    # ------------------------------------------------------------
    # DINÂMICA DOS CAÇADORES (PRODUTORES)
    # ------------------------------------------------------------
      if (is_hunting_day == 1 && time_in_forest > 0) {
        animal_prop <- ifelse((Ia + Sa) > 0, Ia / (Ia + Sa), 0)
        
        # Componentes de risco para caçadores
        injury_risk_hunter <- prob_trans_injury * weight_injury * agg_factor * injury_risk #* (1 - protective_equipment)
        processing_risk_hunter <- prob_trans_processing * weight_processing * processing_risk  #* (1 - protective_equipment)
        eating_risk_hunter <- prob_trans_eating * weight_eating * eating_risk 
        fomite_risk_hunter <- prob_trans_fomite * weight_fomites
        vector_risk_hunter <- prob_trans_vector * bite_risk * encounters_day
        
      # HUNTER DIRECT TRANSMISSION
        hunter_direct_transmission <- time_in_forest * hunting_rate * (
            injury_risk_hunter + 
            processing_risk_hunter + 
            eating_risk_hunter
        ) * animal_prop +
        time_in_forest * cont_rate * fomite_risk_hunter * P
        
        # HUNTER VECTOR-BORNE TRANSMISSION
        #vector_seasonal_factor <- 1 + vector_seasonality * sin(2 * pi * (t - 120) / 365)
        hunter_vector_transmission <- time_in_forest * vector_risk_hunter * animal_prop #* vector_seasonal_factor
        
        # HUNTER SPILLOVER
        hunter_direct_spillover <- (1 - exp(-hunter_direct_transmission)) * Sh
        hunter_vector_spillover <- (1 - exp(-hunter_vector_transmission)) * Sh
        
        hunter_total_spillover_ <- (1 - exp(-hunter_direct_transmission)) + (1 - exp(-hunter_vector_transmission))
        
        dSh <-  -(hunter_direct_spillover + hunter_vector_spillover)
        dIh <-  hunter_direct_spillover + hunter_vector_spillover
      }
      
      # ------------------------------------------------------------
      # DINÂMICA DOS FAMILIARES (CONSUMIDORES)
      # ------------------------------------------------------------
      if (hunting_rate > 0) {
        # Proporção de carne infectada compartilhada
        prop_infected <- ifelse((Ia + Sa) > 0, Ia / (Ia + Sa), 0)
        shared_meat <- sharing_fraction * hunting_rate * prop_infected
        
        # Riscos específicos para consumidores
        processing_risk_consumer <- shared_meat * home_processing_risk #* (1 - home_ppe_effect)
        eating_risk_consumer <- shared_meat * home_eating_risk #* (1 - cooking_protection)
        
        consumer_direct_transmission <- processing_risk_consumer + eating_risk_consumer
        #consumer_spillover_prob <- (1 - exp(-consumer_direct_transmission))
        
        consumer_direct_spillover <- (1 - exp(-consumer_direct_transmission)) * Sc
          
        dSc <-  - consumer_direct_spillover
        dIc <-  consumer_direct_spillover
      }
      
      # ------------------------------------------------------------
      # DINÂMICA ANIMAL E AMBIENTAL
      # ------------------------------------------------------------
      # Transmissão animal-animal
      animal_transmission <- ifelse((Ia + Sa) > 0, 
                                    lambda_effective * Sa * (Ia / (Ia + Sa)), 0)
      
      # Migração
      #migration_in <- migration_rate * K_max * 0.1  # 10% dos migrantes são infectados
      #migration_out <- migration_rate * (A + Ia)
      
      dSa <- r * Sa * (1 - (Sa + Ia) / K_effective) - Sa * (hunting_rate + mu) - 
      animal_transmission #hunting_rate * A - mu * A - 
        #animal_transmission #+ migration_in * 0.9 - migration_out * (A / (A + Ia + 1e-10))
      
      dIa <- animal_transmission - Ia * (hunting_rate + mu) #* Ia - mu * Ia #+ 
        #migration_in * 0.1 - migration_out * (Ia / (A + Ia + 1e-10))
      
      dP <- phi * Ia * (1 - (P / P_max)) - nu * P
      
      # ------------------------------------------------------------
      # OUTPUTS ADICIONAIS PARA ANÁLISE
      # ------------------------------------------------------------
      list(
        c(dSh = dSh, dIh = dIh, dSc = dSc, dIc = dIc, 
          dSa = dSa, dIa = dIa, dP = dP),
        
        # Saídas para análise
        hunter_total_spillover_ = ifelse(exists("hunter_total_spillover_"), hunter_total_spillover_, 0),
        hunter_direct_spillover = ifelse(exists("hunter_direct_spillover"), hunter_direct_spillover, 0),
        hunter_vector_spillover = ifelse(exists("hunter_vector_spillover"), hunter_vector_spillover, 0),
        consumer_direct_spillover = ifelse(exists("consumer_direct_spillover"), consumer_direct_spillover, 0),
        prev_animal_infected = ifelse((Ia + Sa) > 0, Ia / (Ia + Sa), 0),
        hunting_activity = hunting_rate,
        K_effective = K_effective,
        agg_factor = agg_factor,
        lambda_effective = lambda_effective,
        weekday = weekday,
        is_hunting_day = is_hunting_day,
        time_in_forest = time_in_forest,
        climate_effect = climate_effect 
      )
    })
}

#====SIMULACAO====
out <- ode(y = state, times = tempo, func = spillover_model, 
           parms = pars, p_data = precipitation_data$p_normalizado)

# Process results
resultado <- as.data.frame(out) %>%
  left_join(precipitation_data, by = c("time" = "day"))

head(resultado)

#====PLOTS====

# Precipitação
{p1 <- precipitation_data %>% 
  ggplot(aes(x = day, colour = event, group = 1)) +
  geom_point(aes(y = precipitation),
           alpha = 0.6) +
  labs(y = "Daily precipitation (mm)", x = "Day", col = "") +
  theme_bw() + theme(legend.position = "top") +
  scale_color_manual(values = "steelblue")

p2 <- precipitation_data %>%
  ggplot(aes(x = month, y = precipitation)) +
  geom_col(fill = "steelblue") +
  labs(y = "Monthly precipitation (mm)", x = "Month") +
  theme_bw()

precip.p <- ggarrange(
  p1,
  p2,
  ncol = 1,
  common.legend = TRUE
  )
ggsave(
  "./hunt/plots/precip_p.png",
  precip.p,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 300
  )
}

# Variáveis que dependem da precipitação
{
weather_vars <- c(
    "lambda_effective",
    "K_effective",
    "time_in_forest",
    "hunting_activity",
    "climate_effect"
    )

plot_list <- list()

for (var in weather_vars) {
  p <- ggplot(resultado |> filter(is_hunting_day == 1), aes(x = time)) +
    geom_point(aes(y = !!sym(var)),  # Use !!sym() to evaluate string as variable
           color = "steelblue") + 
         labs(title = paste("Time series of", var)) +
    theme_bw()
  

  plot_list[[var]] <- p


weather.p <- ggarrange(
    plot_list$lambda_effective,
    plot_list$K_effective,
    plot_list$time_in_forest,
    plot_list$hunting_activity,
    plot_list$climate_effect,
    ncol = 2,
    nrow = 3
)

ggsave(
  "./hunt/plots/weather_p.png",
  weather.p,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 300
  )
  }
}

# População Animal  
{
animal_pop <- resultado %>%
  ggplot() +
  geom_line(aes(x = time, y = Sa, col = "Susceptible"), linewidth = 0.7) + 
  geom_line(aes(x = time, y = Ia, col = "Infected"), linewidth = 0.7) +
  labs(y = "Animal Population", x = "Time (days)", col = "",
       subtitle = paste("Initial: Sh =", state[1], ", A =", state[2], 
                       ", Ia =", state[3])) +
  scale_color_manual(values = c("Susceptible" = "blue3", "Infected" = "red3")) +
  theme_bw() + theme(legend.position = "top")

ggsave(
  "./hunt/plots/animal_pop.png",
  animal_pop,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 300
  )

}
# População Humana
{

human_pop <- resultado |> select(
  time,Sh,Ih,Sc,Ic) |>
  pivot_longer(cols = 2:5, names_to = "pop", values_to = "n") |> 
  mutate(type = case_when(
    pop %in% c("Ih","Ic") ~ "Infected",
    TRUE ~ "Susceptible"
  ),pop = case_when(
    pop %in% c("Sh","Ih") ~ "Hunters",
    TRUE ~ "Consumers"
  )) |>
  ggplot() +
  geom_line(aes(x = time, y = n, col = type, linetype = pop), linewidth = 0.7) + 
  #geom_line(aes(x = time, y = Ih, col = "Infected"), linewidth = 0.7) +
  labs(y = "Animal Population", x = "Time (days)", col = "", linetype = "",
       subtitle = paste("Initial: Sh =", state[1], ", A =", state[2], 
                       ", Ia =", state[3])) +
  scale_color_manual(values = c("Susceptible" = "blue3", "Infected" = "red3")) +
  theme_bw() + theme(legend.position = "top")

ggsave(
  "./hunt/plots/human_pop.png",
  human_pop,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 300
  )

}
