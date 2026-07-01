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

# NOVA FUNCIONALIDADE: Configuração de cenários
scenario_config <- list(
  baseline = list(name = "Baseline", hunt_intensity = 1.0, vector_season = FALSE),
  high_hunting = list(name = "High Hunting", hunt_intensity = 2.0, vector_season = FALSE),
  vector_season = list(name = "Vector Season", hunt_intensity = 1.0, vector_season = TRUE),
  combined = list(name = "High Hunt + Vectors", hunt_intensity = 2.0, vector_season = TRUE)
)

# Escolher cenário atual
current_scenario <- scenario_config$baseline

# Model parameters - MELHORADO: Parâmetros mais realistas
pars <- list(
  
  # HUNTING
  h = 0.33,#0.14,#0.000001, #0.1 * current_scenario$hunt_intensity,  # hunting efficiency (scenario-dependent)
  h_sazo = 1.2, #0.01,             # seasonal increase 
  base_time_forest = 1,#4,      # average daily hunting time (hours)
  time_max_optimal = 1,#8,      # increased max hunting time
  hunt_days = c(6,7),        # weekend hunting 
  hunter_experience = 0,     # experience factor
  
  
  # ANIMALS
  r = 0.001,#0.001,                  # daily growth rate
  K_max = 50000,#5000,              # maximum carrying capacity
  K_sazo = 0.85,             # seasonal decrease
  mu = 0.0001,                 # daily mortality
  lambda_base = 0.1,        # baseline transmission rate
  lambda_agg_boost = 2.5,    # stronger aggregation effect
  migration_rate = 0,        # migration rate  
  
  # PATHOGENS
  phi = 50,                   # shedding rate
  nu = 0.7,                   # slower pathogen decay
  P_max = 1e6,                # max pathogen load
  pathogen_seasonality = 0.3, #seasonal pathogen survival
  
  # SPILLOVER - Direct transmission (pathogen 1) - for hunters
  injurie_risk = 0, p_I = 0, w_injuries = 0,     #  injury risk
  processing_risk = 0, p_P = 0, w_processing = 0, #  processing risk
  eating_risk = 0, p_E = 0, w_eating = 0,         #  eating risk
  cont_rate = 0, p_F = 0, w_fomites = 0,          #  fomite risk
  protective_equipment = 0, #  PPE usage rate
  
  ## for consumers
  
  sharing_fraction = 0,#0.75,   #0.75,     # meat shared fraction
  home_processing_risk = 0, # home processing risk
  home_eating_risk = 0,#0.1,  #0.25,   # home eating risk
  home_ppe_effect = 0,      # PPE at home
  cooking_protection = 0,   # risk reduction by cooking
  
  
  # SPILLOVER - Vector transmission (pathogen 2)
  encounters_day = 10, #ifelse(current_scenario$vector_season, 10, 10), # seasonal vectors (pensando em usar)
  p_V = 0.15,               #  higher infection probability
  bite_risk = 0.49,          #  higher bite risk
  vector_seasonality = 0    #  seasonal vector activity
)

# Initial conditions
state <- c(
  Sh = 100, #100,               # Hunters S
  Ih = 0,                 # Hunters I
  Sc = 400,               # Consumers S
  Ic = 0,                 # Consumers I
  A = 40000,#4000,                # Animals S
  Ia = 4000,#400,               # Animals I
  P = 0                   # Pathogen in the enviroment
)



spillover_model <- function(t, state, pars, p_data) {
  with(as.list(c(state, pars)), {
    
    # Time indexing
    t_int <- pmax(1, pmin(floor(t), length(p_data)))
    p_norm <- p_data[t_int]
    
    # Climate-dependent effects with more complexity
    threshold <- 0.7  #  lower threshold for more sensitivity
    clim_effect <- pmax(0, 0.3 * (p_norm - threshold) / (1 - threshold))
    agg_factor <- 1 + clim_effect
    
    #  Seasonal pathogen survival
    #pathogen_survival_factor <- 1 + pathogen_seasonality * sin(2 * pi * t / 365)
    
    # carrying capacity
    # seasonal_factor <- 0.8 + 0.4 * sin(2 * pi * (t - 90) / 365)  # vector peak in summer  (not used)
    K_effective <- K_max * (K_sazo + (1 - K_sazo) * (1 - p_norm))
    
    lambda_effective <- lambda_base * (1 + (agg_factor - 1) * lambda_agg_boost) # * pathogen_survival_factor
    
    # Hunting schedule with experience factor
    weekday <- (floor(t) %% 7) + 1  
    is_hunting_day <- ifelse(weekday %in% hunt_days, 1, 0)
    
    #  Weather-dependent hunting (hunters avoid rain)
    weather_factor <- ifelse(p_norm > 0.8, 0, 1.0)  # stop hunting in heavy rain
    
    #  More complex time in forest calculation
    optimal_climate <- 0.2  # optimal precipitation level for hunting
    climate_hunting_factor <- 1 - 0.5 * abs(p_norm - optimal_climate)
    
    time_in_forest_base <- (base_time_forest + 
                              (time_max_optimal - base_time_forest) * 
                              climate_hunting_factor) / 24 #* hunter_experience) / 24
    
    # Effective rates on hunting days only
    time_in_forest <- is_hunting_day * time_in_forest_base * weather_factor
    hunting_rate <- is_hunting_day * h * (1 + h_sazo * p_norm) #* hunter_experience
    
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
      ) * animal_prop #+ fomite_risk
      
      # Transmissão vetorial
      vector_seasonal_factor <- 1 + vector_seasonality * sin(2 * pi * (t - 120) / 365)
      vector_risk_hunter <- time_in_forest * encounters_day * 
        bite_risk * p_V * animal_prop * vector_seasonal_factor
      
      hunter_spillover <- (1 - exp(-hunter_risk)) * Sh #* (Ia/A+Ia)
      vector_spillover <- (1 - exp(-vector_risk_hunter)) * Sh
      
      hunter_spillover_prob <- (1 - exp(-hunter_risk))+ (1 - exp(-vector_risk_hunter))
      
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
      
      consumer_spillover_prob <- (1 - exp(-consumer_risk))
      
      consumer_spillover <- (1 - exp(-consumer_risk)) * Sc
      
      dSc <-  - consumer_spillover
      dIc <-  consumer_spillover
    }
    
    # ------------------------------------------------------------
    # DINÂMICA ANIMAL E AMBIENTAL
    # ------------------------------------------------------------
    # Transmissão animal-animal
    animal_transmission <- ifelse((Ia + A) > 0, 
                                  lambda_effective * A * (Ia / (Ia + A)), 0)
    
    # Migração
    migration_in <- migration_rate * K_max * 0.1  # 10% dos migrantes são infectados
    migration_out <- migration_rate * (A + Ia)
    
    # Equações diferenciais
    dA <- r * A * (1 - (A + Ia) / K_effective) - Sh * hunting_rate * A/(Ia+A) - mu * A - 
      animal_transmission +
      #migration_in * 0.9 - migration_out * (A / (A + Ia + 1e-10)) + 
      0.1 * Ia
    
    dIa <- animal_transmission - Sh * hunting_rate * Ia/(Ia+A) - (mu*5) * Ia - 
      #migration_in * 0.1 - migration_out * (Ia / (A + Ia + 1e-10)) -
      0.1 * Ia
    
    dP <- phi * Ia * (1 - (P / P_max)) - nu * P
    
    # ------------------------------------------------------------
    # OUTPUTS ADICIONAIS PARA ANÁLISE
    # ------------------------------------------------------------
    list(
      c(dSh = dSh, dIh = dIh, dSc = dSc, dIc = dIc, 
        dA = dA, dIa = dIa, dP = dP),
      
      # Saídas para análise
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
      migration_in = migration_in,
      migration_out = migration_out, 
      hunter_spillover_prob = hunter_spillover_prob,
      consumer_spillover_prob = consumer_spillover_prob
      
    )
  })
}

# Run simulation
out <- ode(y = state, times = tempo, func = spillover_model, 
           parms = pars, p_data = precipitation_data$p_normalizado)

# Process results
resultado <- as.data.frame(out) %>%
  left_join(precipitation_data, by = c("time" = "day")) #%>%
  # mutate(Sh = (Sh/100) * 100,
  #        prev = Ia / (A + Ia))

resultado %>%
  ggplot() +
  geom_line(aes(x = time, y = A, col = "Susceptible"), linewidth = 0.7) +
  geom_line(aes(x = time, y = Ia, col = "Infected"), linewidth = 0.7) +
  labs(y = "Animal Population", x = "Time (days)", col = "",
       subtitle = paste("Initial: Sh =", state[1], ", A =", state[2],
                        ", Ia =", state[3])) +
  scale_color_manual(values = c("Susceptible" = "blue3", "Infected" = "red3")) +
  theme_bw() + theme(legend.position = "top")
# #
# #
# resultado %>%
#   ggplot()+
#   geom_line(aes(x=time, y=prev),col="orange3")

# 
# resultado.sim1 <- resultado %>% select(time,Sh,Ih,Sc,Ic,A,Ia) %>%
#   mutate(animal = "D. rotundus")

# resultado.sim2 <- resultado %>% select(time,Sh,Ih,Sc,Ic,A,Ia) %>%
#      mutate(animal = "A. seniculus")

resultado.sim3 <- resultado %>% select(time,Sh,Ih,Sc,Ic,A,Ia) %>%
       mutate(animal = "C. paca")
# 
# resultado.simAnimal4 <- rbind(resultado.sim1,resultado.sim2,resultado.sim3)
# #
#  resultado.simAnimal4 <- resultado.simAnimal4 %>%
#    mutate(prev = Ia/(A+Ia),
#           year = ((time-1)%/%365)+1)
#
# prevIa <- ggplot(resultado.simAnimal) +
#   geom_line(aes(x=time,y=prev*100, col = animal), alpha = 0.8) +
#   theme_bw(base_size = 12) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Prevalence of infected animals (%)") +
#   theme(legend.text = element_text(face = "italic")) ; prevIa

# ggsave("./book_chapter/prevIa2.5.png",prevIa,dpi = 300, width = 6,height = 4)
# 
# dynaA_Ia_hunt <- ggplot(resultado.simAnimal) +
#   geom_line(aes(x=time,y=A, col = animal)) +
#   geom_line(aes(x=time,y=Ia, col = animal)) +
#   theme_bw(base_size = 12) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Number of susceptible animals") +
#   theme(legend.text = element_text(face = "italic")); dynaA_Ia_hunt

# ggsave("./book_chapter/dynaA_Ia_hunt.png",dynaA_hunt,dpi = 300, width = 6,height = 4)


# ggplot(resultado.simAnimal) +
#   geom_line(aes(x=time,y=A, col = animal)) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#     scale_x_continuous(
#       breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#       labels = c(seq(1,10,1))
#     ) +
#     labs(x = "Time (years)", y = "Number of susceptible animals") +
#     theme(legend.text = element_text(face = "italic"))
# 
#  
#  sim1.figA <-  ggplot(resultado.simAnimal) +
#    geom_line(aes(x=time,y=A, col = animal)) +
#    #geom_line(aes(x=time,y=Ia, col = animal), linetype = 2) +
#    scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#    scale_x_continuous(
#      breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#      labels = c(seq(1,10,1))
#    ) +
#    labs(x = "Time (years)", y = "Number of animals") +
#    theme_bw(base_size = 12) +
#    theme(legend.text = element_text(face = "italic")); sim1.figA
# 
# ggsave("./book_chapter/animal_cenarioA.png",sim1.figA, dpi = 300, bg = "white",
#        width = 6,height = 4)
 
 
# # 
# sim1.figB <-  ggplot(resultado.simAnimal2) +
#   geom_line(aes(x=time,y=A, col = animal)) +
#   geom_line(aes(x=time,y=Ia, col = animal), linetype = 2, linewidth = 0.5) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Number of animals") +
#   theme_bw(base_size = 12) +
#   theme(legend.text = element_text(face = "italic")) ; sim1.figB
# 
# ggsave("./book_chapter/animal_cenarioB.png",sim1.figB, dpi = 300, bg = "white",
#               width = 6,height = 4)

# 
# sim1.figC <-  ggplot(resultado.simAnimal3) +
#   geom_line(aes(x=time,y=A, col = animal)) +
#   #geom_line(aes(x=time,y=Ia, col = animal), linetype = 2) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Number of animals") +
#   theme_bw(base_size = 12) +
#   theme(legend.text = element_text(face = "italic")); sim1.figC
# 
# ggsave("./book_chapter/animal_cenarioC.png",sim1.figC, dpi = 300, bg = "white",
#                      width = 6,height = 4)

# 
# sim1.figD <-  ggplot(resultado.simAnimal4) +
#   geom_line(aes(x=time,y=A, col = animal)) +
#   geom_line(aes(x=time,y=Ia, col = animal), linetype = 2, linewidth = 0.5) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Number of animals") +
#   theme_bw(base_size = 12) +
#   theme(legend.text = element_text(face = "italic")) ; sim1.figD
# 
# ggsave("./book_chapter/animal_cenarioD.png",sim1.figD, dpi = 300, bg = "white",
#                      width = 6,height = 4)

# 
# sim1.figAB  <- ggarrange(sim1.figA + labs(x="",y = ""),
#           sim1.figB + labs(x="",y = ""),
#           # sim1.figC + labs(x="",y = ""),
#           # sim1.figD + labs(x="",y = ""),
#           labels = letters[1:4],align = "v",
#           ncol = 1, nrow = 2, common.legend = T) %>%
#   annotate_figure(
#                   left = text_grob("Number of animals",
#                                    rot = 90, size = 14),
#                   bottom = text_grob("Time (years)",
#                                      size = 14),); sim1.figAB
# 
# ggsave("./book_chapter/sim1.figAB.png",sim1.figAB,dpi = 300, width = 6,height = 7, bg = "white")

# sim1.figCD  <- ggarrange(sim1.figC + labs(x="",y = ""),
#           sim1.figD + labs(x="",y = ""),
#           # sim1.figC + labs(x="",y = ""),
#           # sim1.figD + labs(x="",y = ""),
#           labels = letters[1:4],align = "v",
#           ncol = 1, nrow = 2, common.legend = T) %>%
#   annotate_figure(
#                   left = text_grob("Number of animals",
#                                    rot = 90, size = 14),
#                   bottom = text_grob("Time (years)",
#                                      size = 14),); sim1.figCD

# ggsave("./book_chapter/sim1.figCD.png",sim1.figCD,dpi = 300, width = 6,height = 7, bg = "white")

# sim1.figABCD  <- ggarrange(sim1.figA + labs(x="",y = ""),
#                          sim1.figB + labs(x="",y = ""),
#                          sim1.figC + labs(x="",y = ""),
#                          sim1.figD + labs(x="",y = ""),
#                          labels = letters[1:4], align = "v",
#                          ncol = 2, nrow = 2, common.legend = T) %>%
#   annotate_figure(
#     left = text_grob("Number of animals",
#                      rot = 90, size = 14),
#     bottom = text_grob("Time (years)",
#                        size = 14),); sim1.figABCD  
# 
# ggsave("./book_chapter/sim1.figABCD.png",sim1.figABCD,dpi = 300, width = 9,height = 7, bg = "white")


# 
# 
# sim2.figB <- ggplot(resultado.simAnimal2) +
#   geom_line(aes(x=time,y=prev*100, col = animal), alpha = 1, linetype = 1) +
#   theme_bw(base_size = 12) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Prevalence of infected animals (%)") +
#   theme(legend.text = element_text(face = "italic")) ; sim2.figB
# 
# sim2.figD <- ggplot(resultado.simAnimal4) +
#   geom_line(aes(x=time,y=prev*100, col = animal), alpha = 1, linetype = 1) +
#   theme_bw(base_size = 12) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Prevalence of infected animals (%)") +
#   theme(legend.text = element_text(face = "italic")) ; sim2.figD

sim2.figBD  <- ggarrange(sim2.figB + labs(x="",y = "") + ylim (0,15),
                         sim2.figD + labs(x="",y = "") + ylim (0,15),
          # sim1.figC + labs(x="",y = ""),
          # sim1.figD + labs(x="",y = ""),
          labels = letters[1:4],align = "v",
          ncol = 1, nrow = 2, common.legend = T) %>%
  annotate_figure(
                  left = text_grob("Prevalence (%)",
                                   rot = 90, size = 14),
                  bottom = text_grob("Time (years)",
                                     size = 14),); sim2.figBD
#  
  ggsave("./book_chapter/sim2.figBD.png",sim2.figBD,dpi = 300, width = 6,height = 7, bg = "white")

# 
# 
# sim2.figB <- ggplot(resultado.simAnimal4) +
#   geom_line(aes(x=time,y=prev*100, col = animal), alpha = 1, linetype = 2) +
#   theme_bw(base_size = 12) +
#   scale_color_manual(values = c("forestgreen","chocolate3","firebrick3")) +
#   scale_x_continuous(
#     breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#     labels = c(seq(1,10,1))
#   ) +
#   labs(x = "Time (years)", y = "Prevalence of infected animals (%)") +
#   theme(legend.text = element_text(face = "italic")) ; sim2.figB
# 
# sim2.fig  <- ggarrange(sim2.figA + labs(x="",y = ""),
#           sim2.figB + labs(x="",y = ""),
#           labels = letters[1:2],
#           ncol = 1, nrow = 2, common.legend = T) %>%
#   annotate_figure(
#                   left = text_grob("Prevalence of infected animals (%)",
#                                    rot = 90, size = 14),
#                   bottom = text_grob("Time (years)",
#                                      size = 14)) ; sim2.fig
# 
# ggsave("./book_chapter/sim2.fig.png",sim2.fig,dpi = 300, width = 9,height = 6, bg = "white")

#### simulation with hunting transmission

# scenario C (no infected but hunting)

# resultado.sim1 %>% 
#   mutate(prev_hum = Ih/(Sh+Ih),
#          prev_anim = Ia/(Ia+A)) %>% 
#   ggplot() +
#   geom_line(aes(x = time, y = prev_hum)) + 
#   geom_line(aes(x = time, y = prev_anim))
# 
# # scenario D (infected with hunting)
# 

# sim3.rabies <- resultado %>%
#   mutate(prev_hun = Ih/(Sh+Ih),
#          prev_anim = Ia/(Ia+A),
#          prev_con = ifelse(is.na(Ic/(Ic+Sc)),0,Ic/(Ic+Sc))) %>%
#   ggplot() +
#   geom_line(aes(x = time, y = prev_anim *100, col = "animal"), linetype = 1) +
#   geom_line(aes(x = time, y = prev_hun * 100, col = "hunter"), linetype = 1) +
#   geom_line(aes(x = time, y = prev_con * 100, col = "consumers"), linetype = 1) +
#   scale_color_manual(values = c("red","#1E90FF","black")) +
#   theme_bw(base_size = 15) + 
#   theme(legend.position = "top") +
#   labs(x = "Time (years)", y = "Prevalence (%)", color = "") +
#   scale_x_continuous(
#         breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#         labels = c(seq(1,10,1))
#       ); sim3.rabies 
# 
# ggsave("./book_chapter/sim3.rabies.png",sim3.rabies,dpi = 300, width = 6,height = 7, bg = "white")

# sim3.vogeli <- resultado %>%
#     mutate(prev_hun = Ih/(Sh+Ih),
#            prev_anim = Ia/(Ia+A),
#            prev_con = ifelse(is.na(Ic/(Ic+Sc)),0,Ic/(Ic+Sc))) %>%
#     ggplot() +
#     geom_line(aes(x = time, y = prev_anim *100, col = "animal"), linetype = 1) +
#     geom_line(aes(x = time, y = prev_hun * 100, col = "hunter"), linetype = 1) +
#     geom_line(aes(x = time, y = prev_con * 100, col = "consumers"), linetype = 1) +
#     scale_color_manual(values = c("red","#1E90FF","black")) +
#     theme_bw(base_size = 15) +
#     theme(legend.position = "top") +
#     labs(x = "Time (years)", y = "Prevalence (%)", color = "") +
#     scale_x_continuous(
#           breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
#           labels = c(seq(1,10,1))
#         ); sim3.vogeli
# 
# ggsave("./book_chapter/sim3.vogeli.png",sim3.vogeli,dpi = 300, width = 6,height = 7, bg = "white")
# 

sim3.mayaro <- resultado %>%
    mutate(prev_hun = Ih/(Sh+Ih),
           prev_anim = Ia/(Ia+A),
           prev_con = ifelse(is.na(Ic/(Ic+Sc)),0,Ic/(Ic+Sc))) %>%
    ggplot() +
    geom_line(aes(x = time, y = prev_anim *100, col = "animal"), linetype = 1) +
    geom_line(aes(x = time, y = prev_hun * 100, col = "hunter"), linetype = 1) +
    geom_line(aes(x = time, y = prev_con * 100, col = "consumers"), linetype = 1) +
    scale_color_manual(values = c("red","#1E90FF","black")) +
    theme_bw(base_size = 15) +
    theme(legend.position = "top") +
    labs(x = "Time (years)", y = "Prevalence (%)", color = "") +
    scale_x_continuous(
          breaks = seq(1, max(resultado.simAnimal$time), by = 365),  # Breaks at 0, 365, 730
          labels = c(seq(1,10,1))
        ); sim3.mayaro

ggsave("./book_chapter/sim3.mayaro.png",sim3.mayaro,dpi = 300, width = 6,height = 7, bg = "white")


sim3.figRVM <- ggarrange(sim3.rabies + labs(x="",y = "") + ylim(0,100),
                            sim3.vogeli + labs(x="",y = "") + ylim(0,100),
                            sim3.mayaro + labs(x="",y = "") + ylim(0,100),
                            labels = letters[1:4], align = "v",
                            ncol = 3, nrow = 1, common.legend = T) %>%
                             annotate_figure(
                               left = text_grob("Prevalence (%)",
                                                rot = 90, size = 14),
                               bottom = text_grob("Time (years)",
                                                  size = 14),); sim3.figRVM

ggsave("./book_chapter/sim3.figRVM.png",sim3.figRVM,dpi = 300, width = 9,height = 7, bg = "white")

