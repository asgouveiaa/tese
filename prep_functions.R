## Functions for hunting model

generate_precipitation_dataset <- function(
    years = NULL,                     # Número de anos a simular
    days_per_year = 365,
    base = 1,                      # Precipitação base
    amplitude = 3,                 # Amplitude da variação sazonal
    noise_sd = 0.5,                # Intensidade do ruído
    la_nina_years = NULL,          # Anos com La Niña (vetor numérico)
    nina_intensity = 0,            # Intensidade do efeito La Niña
    drought_periods = NULL,        # Períodos de seca (lista com start/end)
    flood_periods = NULL,          # Períodos de enchente (lista com start/end)
    flood_intensity = 5,           # Fator de intensidade de enchentes
    random_seed = NULL             # Semente para reproducibilidade
) {
  
  # Configurar semente aleatória se fornecida
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  # Total de dias a simular
  total_days <- years * 365
  
  # Inicializar dataframe
  precipitation_data <- data.frame(
    day = 1:total_days,
    year = ceiling((1:total_days) / 365),
    year_day = ((1:total_days) - 1) %% 365 + 1,
    precipitation = numeric(total_days),
    is_la_nina = numeric(total_days),
    is_drought = numeric(total_days),
    is_flood = numeric(total_days),
    stringsAsFactors = FALSE
  )
  
  # Gerar dados para cada dia
  for (i in 1:total_days) {
    # Calcular efeitos sazonais
    phase <- -pi / 2
    seasonal_precip <- base + amplitude * (1 + sin(2 * pi * precipitation_data$year_day[i] / 365 + phase))
    
    # Verificar La Niña
    is_la_nina <- ifelse(precipitation_data$year[i] %in% la_nina_years, 1, 0)
    
    # Aplicar efeito La Niña
    if (is_la_nina == 1) {
      nina_effect <- ifelse(precipitation_data$year_day[i] >= 274 | precipitation_data$year_day[i] <= 90, 
                            1 + nina_intensity, 1)
      seasonal_precip <- seasonal_precip * nina_effect
    }
    
    # Verificar seca
    is_drought <- FALSE
    if (!is.null(drought_periods)) {
      for (period in drought_periods) {
        if (i >= period$start && i <= period$end) {
          seasonal_precip <- 0
          is_drought <- TRUE
          break
        }
      }
    }
    
    # Verificar enchente
    is_flood <- FALSE
    if (!is_drought && !is.null(flood_periods)) {
      for (period in flood_periods) {
        if (i >= period$start && i <= period$end) {
          seasonal_precip <- seasonal_precip * flood_intensity
          is_flood <- TRUE
          break
        }
      }
    }
    
    # Adicionar ruído (exceto durante secas)
    if (!is_drought) {
      precipitation <- seasonal_precip + rnorm(1, mean = 0, sd = noise_sd)
    } else {
      precipitation <- 0
    }
    
    # Armazenar resultados
    precipitation_data$precipitation[i] <- max(0, precipitation)
    precipitation_data$is_la_nina[i] <- is_la_nina
    precipitation_data$is_drought[i] <- as.numeric(is_drought)
    precipitation_data$is_flood[i] <- as.numeric(is_flood)
  }
  
  # Adicionar coluna de data real (assumindo início em 1º de janeiro de 2000)
  precipitation_data$date <- as.Date("2000-01-01") + (precipitation_data$day - 1)
  
  # Adicionar coluna de mês
  precipitation_data$month <- format(precipitation_data$date, "%m")
  
  # Adicionar coluna tipo de evento
  precipitation_data %>% 
    mutate(event = case_when(
      is_drought == 1 ~ "drought",
      is_flood == 1 ~ "flood",
      TRUE ~ "normal"
    )) -> precipitation_data
  
  # Adicionar coluna La_nina
  precipitation_data %>% 
    mutate(la_nina_year = case_when(
      is_la_nina == 1 ~ "yes",
      TRUE ~ "no"
    )) -> precipitation_data
  
  # Dataset final
  precipitation_data %>% 
    select(day,year,year_day,date,month,event,
           la_nina_year,precipitation) -> precipitation_data
  
  return(precipitation_data)
}

# Função auxiliar para criar períodos mais facilmente
create_periods <- function(start_days, end_days) {
  if (length(start_days) != length(end_days)) {
    stop("start_days e end_days devem ter o mesmo comprimento")
  }
  periods <- list()
  for (i in seq_along(start_days)) {
    periods[[i]] <- list(start = start_days[i], end = end_days[i])
  }
  return(periods)
}


# calculate_drainage <- function(precipitation_daily, lambda = lambda, threshold = threshold) {
#   n <- length(precipitation_daily)
#   funS <- numeric(n)  # Vetor para armazenar a saturação diária
  
#   for (t in 1:n) {
#     # Soma acumulada com decaimento exponencial
#     decay_effect <- sapply(1:t, function(i) precipitation_daily[i] * exp(-lambda * (t - i)))
#     funS[t] <- min(threshold, sum(decay_effect))
#   }
  
#   return(funS)
# }