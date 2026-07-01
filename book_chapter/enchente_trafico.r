ana <- read.csv2("~/Documents/git/spillover25/ana_fonteboa.csv")

head(ana)
glimpse(ana)

ana %>% mutate(mes = substr(Data,4,5),
               ano = substr(Data,7,10)) %>% 
  filter(hora == "07:00") %>% 
  select(ano,mes,starts_with("Cota")) %>% 
  pivot_longer(cols = starts_with("Cota"),names_to = "day", values_to = "mm") %>% 
  mutate(day = as.numeric(substr(day,5,6))) -> ana2

ana2 <- ana2 %>% 
  mutate(
    data = make_date(ano, mes, day),
    dia_do_ano = yday(data),
    target = ifelse(ano <= 202,0,1)
  ) 

ana2 %>% 
  ggplot(aes(dia_do_ano, mm/100, color = factor(target), group = ano)) +
  geom_line()

ana2 %>% 
  mutate(ano = as.numeric(ano)) %>% 
  filter(ano %in% c(1979:2022)) %>% 
  select(ano,dia_do_ano,mm) %>% 
  group_by(dia_do_ano) %>% 
  summarise(mean = mean(mm,na.rm = T),
            sd = sd(mm)) %>% 
  mutate(lim_s = mean + 1.96*(sd/43),
         lim_i = mean - 1.96*(sd/43)) -> ana_resumo


# Filtrar período de referência (1979-2022)
ana_referencia <- ana2 %>%
  filter(ano >= 1979 & ano <= 2022)

# Calcular estatísticas por dia do ano
ana_resumo <- ana2 %>%
  select(dia_do_ano,ano,mm) %>% 
  group_by(dia_do_ano) %>%
  summarise(
    media = mean(mm, na.rm = TRUE),
    desvio_padrao = sd(mm, na.rm = TRUE),
    n = n(),
    erro_padrao = desvio_padrao / sqrt(n),
    ic_inferior = media - (1.96 * erro_padrao),
    ic_superior = media + (1.96 * erro_padrao),
    .groups = 'drop'
  )

ana_resumo %>% 
  ggplot()+
  geom_line(aes(dia_do_ano,media/100)) +
  geom_ribbon(aes(dia_do_ano,ymin = ic_inferior/100,
                  ymax = ic_superior/100), alpha = 0.3)



# data plotly via deepseek ------------------------------------------------

ana <- read.csv2("~/Documents/git/spillover25/enso_data_completo.csv")
head(ana)

glimpse(ana)

enchente.fig <- ana %>% 
  mutate(dia = 1:nrow(ana)) %>% 
  ggplot() + 
  # Mapeia a linha tracejada na estética de cor
  geom_line(aes(dia, as.numeric(Mean.Média.1979.2022), linetype = "Mean 1979-2022"), size = 0.8) +
  # Mapeia a ribbon na estética de preenchimento
  geom_ribbon(aes(dia, ymin = as.numeric(X95..CI...IC.1),
                  ymax = as.numeric(X95..CI...IC), fill = "95% CI"), alpha = 0.3) + 
  geom_line(aes(dia, as.numeric(X2022), color = "2022")) +
  geom_line(aes(dia, as.numeric(X2023), color = "2023")) +
  geom_line(aes(dia, as.numeric(X2024), color = "2024")) +
  geom_line(aes(dia, as.numeric(X2025), color = "2025")) + 
  scale_color_manual(
    name = "Years",
    values = c(
      "2022" = "#1f77b4",  # Azul
      "2023" = "#ff7f0e",  # Laranja
      "2024" = "#2ca02c",  # Verde
      "2025" = "#d62728"   # Vermelho
    )
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Mean 1979-2022" = 2)
  ) +
  scale_fill_manual(
    name = "",
    values = c("95% CI" = "gray50")
  ) +
  labs(
    title = "",
    y = "Water level (m)", 
    x = "Day of the year"
  ) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top") #+ 
  #scale_x_discrete(brseq(0,365,1),expand = c(0,0)); enchente.fig


ggsave("./figs/enchente.fig.png", enchente.fig ,dpi = 600, 
       width = 9, height = 5, bg = "white")


## tráfico de animais

data <- data.frame(
  group = c("Reptiles", "Reptiles", "Reptiles", "Reptiles", "Reptiles", 
            "Birds", "Birds", "Birds", "Birds", "Birds",
            "Mammals", "Mammals", "Mammals", "Mammals", "Mammals",
            "Fish", "Fish", "Fish", "Fish", "Fish",
            "Insects", "Insects", "Insects", "Insects", "Insects",
            "Amphibians", "Amphibians", "Amphibians", "Amphibians", "Amphibians",
            "Arthropods", "Arthropods", "Arthropods", "Arthropods", "Arthropods"),
  common_name = c("Amazon river turtle", "Yellow-spotted river turtle", "Iguana", "Boa constrictor", "Rainbow boa",
                  "Thick-billed saltator", "Amazonian canary", "Large-billed seed finch", "Orange-winged amazon", "Blue-and-yellow macaw",
                  "Red howler monkey", "Black-striped capuchin", "Black spider monkey", "Brown-throated sloth", "Amazonian manatee",
                  "Arapaima", "Zebrafish", "Plecostomus", "Pearl ray", "Motor ray",
                  "Longhorn beetle", "Hercules beetle", "Titan beetle", "Blue morpho butterfly", "Owl butterfly",
                  "Amazon frogs", "Poison dart frog", "Horned frog", "Giant monkey frog", "Salamander",
                  "Goliath birdeater", "Brazilian crab spider", "Giant tarantula", "Black scorpion", "Giant centipede"),
  purpose = c("Food", "Food", "Pet", "Pet", "Pet",
              "Pet", "Pet", "Pet", "Traditional use / Pet", "Traditional use / Pet",
              "Food/Pet", "Food/Pet", "Food/Pet", "Pet", "Food",
              "Food", "Hobby (collecting)", "Hobby (collecting)", "Hobby (collecting)", "Hobby (collecting)",
              "Pet", "Pet", "Pet", "Hobby (collecting)", "Hobby (collecting)",
              "Pet", "Pet", "Pet", "Pet", "Pet",
              "Pet", "Pet", "Pet", "Pet", "Pet")
)

data

# Certifique-se de que o pacote ggplot2 está instalado
# install.packages("ggplot2")

dados_percentual <- data %>%
  # Conta o número de espécies para cada GRUPO e FINALIDADE
  group_by(group, purpose) %>%
  summarise(n = n(), .groups = 'drop') %>%
  # Calcula a porcentagem de 'n' em relação ao total de espécies em cada 'GRUPO'
  group_by(group) %>%
  mutate(percentual = n / sum(n))

escala_qualitativa <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")
escala_pastel <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")


# Geração do Gráfico de Barras Empilhadas por Porcentagem
fig.traffic <- ggplot(dados_percentual, aes(x = group, y = percentual*100, fill = purpose)) +
  
  # 1. Barras empilhadas
  geom_bar(stat = "identity", position = "stack", color = "black") +
  
  # 2. Rótulos de porcentagem dentro das barras
  geom_text(
    aes(label = scales::percent(percentual)), # Formata como porcentagem
    position = position_stack(vjust = 0.5), # Centraliza no meio da pilha
    color = "black",
    size = 3.5,
    fontface = "bold"
  ) +
  
  # 3. Escalas e Títulos
  # scale_y_continuous(labels = scales::percent) + # Eixo Y em formato de porcentagem
  # labs(
  #   title = "",
  #   x = "Grupo Taxonômico",
  #   y = "Percentual de Espécies (%)",
  #   fill = "Finalidade"
  # ) +
  # scale_fill_manual(values = escala_pastel) +
  scale_fill_brewer(palette = "Spectral") +  # Para 5 cores automaticamente
  labs(x = "", y = "Percentual of animals (%)", fill = "Purpose") +
  # 4. Tema e Ajustes Visuais
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotaciona os rótulos do eixo X
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  theme(legend.position = "top"); fig.traffic

ggsave("./figs/fig.traffic.png", fig.traffic ,dpi = 1000, 
       width = 7, height = 5, bg = "white")


# 1. Criação do Data Frame
dados_trafico <- data.frame(
  POSICAO = c("1º", "2º", "3º", "4º", "5º", "6º", "7º", "8º", "9º", "10º", "11º", "12º", "13º", "14º", "15º", "16º", "17º", "18º", "19º", "20º"),
  MUNICIPIO = c("BARCELOS", "SANTA ISABEL DO RIO NEGRO", "ALTAMIRA", "MEDICILÂNDIA", "ILHA DO MARAJÓ", "ALMEIRIM", "SERRA DO NAVIO", "MANAUS", "LABREA", "ITACOATIARA", "SANTARÉM", "GURUPI", "ALVORADA", "BOA VISTA", "CRUZEIRO DO SUL", "CARACARAÍ", "PRESIDENTE FIGUEIREDO", "PORTO VELHO", "ARAGUAÍNA", "COLNIZA"),
  UF = c("AM", "AM", "PA", "PA", "PA", "PA", "AP", "AM", "AM", "AM", "PA", "TO", "TO", "RR", "AC", "RR", "AM", "RO", "TO", "MT"),
  ESPECIES = c(
    "PEIXES ORNAMENTAIS", "PEIXES ORNAMENTAIS", "PEIXES ORNAMENTAIS / AVES / RÉPTEIS", "PEIXES ORNAMENTAIS / RÉPTEIS",
    "RÉPTEIS / PEIXES ORNAMENTAIS", "PEIXES ORNAMENTAIS / AVES / FELINOS", "ANFIBIOS / AVES / FELINOS / PRIMATAS",
    "PRIMATAS / FELINOS / AVES", "AVES / PRIMATAS / RÉPTEIS", "QUELÔNIOS / ARACNÍDEOS", "PEIXES ORNAMENTAIS / ARACNÍDEOS",
    "AVES", "AVES / RÉPTEIS", "AVES", "AVES /", "AVES / RÉPTEIS", "AVES / PRIMATAS", "PRIMATAS / RÉPTEIS", "AVES/RÉPTEIS", "FELINOS"
  )
)

# 2. Limpeza e "Tidying" (Transformação de dados)
dados_tidy <- dados_trafico %>%
  # Limpa o separador (para lidar com '/ ' ou '/') e separa as strings em múltiplas linhas
  mutate(ESPECIES = gsub(" / |/", ",", ESPECIES)) %>% 
  separate_rows(ESPECIES, sep = ",") %>%
  # Limpa espaços em branco e remove linhas vazias que podem ter sido geradas
  mutate(GRUPO_TRAFICADO = trimws(ESPECIES)) %>%
  filter(GRUPO_TRAFICADO != "") %>%
  
  # Cria uma coluna que combina Município e UF para o eixo X, e reordena pelo ranking (POSICAO)
  mutate(
    MUNICIPIO_UF = paste(MUNICIPIO, " (", UF, ")", sep=""),
    # Converte POSICAO em número para ordenação correta, depois em factor para o eixo X
    POSICAO_NUM = as.numeric(gsub("[º]", "", POSICAO))
  ) %>%
  arrange(POSICAO_NUM) %>%
  # Define a ordem final do eixo X
  mutate(MUNICIPIO_UF = factor(MUNICIPIO_UF, levels = unique(MUNICIPIO_UF)))

ggplot(dados_tidy, aes(x = MUNICIPIO_UF, fill = GRUPO_TRAFICADO)) +
  
  # 1. Cria o gráfico de barras (as barras são empilhadas por padrão com 'fill')
  geom_bar(color = "black") +
  
  # 2. Rótulos e Títulos
  labs(
    title = "Grupos de Espécies Traficadas por Município (Top 20)",
    subtitle = "Altura da barra indica a diversidade de grupos traficados no local.",
    x = "Município (UF)",
    y = "Contagem de Grupos de Espécies Traficadas",
    fill = "Grupo de Espécies"
  ) +
  
  # 3. Tema e Ajustes Visuais
  scale_fill_brewer(palette = "Spectral") + # Usa uma paleta de cores para clareza
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9), # Rotaciona os rótulos do eixo X
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + coord_flip()

require(gt)

tabela_trafico_gt <- dados_trafico %>%
  gt() %>%
  
  # Adiciona título e subtítulo
  tab_header(
    title = md("**Top 20 Municípios com Maior Registro de Tráfico de Espécies**"),
    #subtitle = "Dados de Grupos de Espécies Traficadas"
  ) %>%
  
  # Renomeia as colunas
  cols_label(
    POSICAO = "POSIÇÃO",
    MUNICIPIO = "MUNICÍPIO",
    UF = "UF",
    ESPECIES = "GRUPOS DE ESPÉCIES"
  ) %>%
  
  # Adiciona linhas de separação no corpo da tabela
  tab_options(
    row.striping.include_table_body = TRUE,
    row.striping.background_color = "#f7f7f7"
  ) %>%
  
  # Aplica um estilo de fonte/tema (opcional, mas recomendado)
  opt_table_font(font = list(c("Arial", "sans-serif"))) %>%
  
  # Alinha o texto das colunas
  cols_align(align = "center", columns = c(POSICAO, UF)) %>%
  cols_align(align = "left", columns = c(MUNICIPIO, ESPECIES)) %>%
  
  # Formata o cabeçalho (posição)
  tab_spanner(
    label = "RANKING",
    columns = POSICAO
  )

tabela_trafico_gt

# 3. Salva a tabela como imagem (PNG)
# O arquivo será salvo na sua pasta de trabalho (working directory)
gtsave(tabela_trafico_gt, filename = "tabela_trafico.png")
