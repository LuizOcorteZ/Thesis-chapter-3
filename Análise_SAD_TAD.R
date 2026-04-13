# ==============================================================================
# SCRIPT DOUTORADO - CAPÍTULO 3 (VERSÃO FINAL CONSOLIDADA)
# OBJETIVO: Respostas de TAD e SAD ao Gradiente de Paisagem e Solo
# ==============================================================================

# 1. CARREGAR PACOTES ----------------------------------------------------------
pacotes <- c("tidyverse", "vegan", "e1071", "ade4", "openxlsx", 
             "broom", "car", "gridExtra", "grid", "cowplot")

for(p in pacotes) {
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# 2. CARREGAR E ALINHAR DADOS --------------------------------------------------
master <- read.xlsx("dados_completos2.xlsx")
abund  <- read.xlsx("abundance_final_itv.xlsx", rowNames = TRUE)
traits <- read.xlsx("traits_final_corrected.xlsx", rowNames = TRUE)

# Limpeza de nomes e match de espécies
colnames(abund) <- gsub("[ _]+", ".", trimws(colnames(abund)))
colnames(abund) <- gsub("\\.\\.+", ".", colnames(abund))
rownames(traits) <- gsub("[ _]+", ".", trimws(rownames(traits)))
rownames(traits) <- gsub("\\.\\.+", ".", rownames(traits))

comuns <- intersect(colnames(abund), rownames(traits))
abund_final  <- abund[, comuns] %>% .[, sort(colnames(.))]
traits_final <- traits[comuns, ] %>% .[order(rownames(.)), ]

stopifnot(all(colnames(abund_final) == rownames(traits_final)))

if(!"uso" %in% colnames(master)) master$uso <- substr(master$cod, 1, 1)



# --- DEBUG DO LOOP ---
i <- 1 # Vamos testar apenas o primeiro sítio (C1a)

print(paste("Analisando sítio:", rownames(abund_final)[i]))

# 1. O que tem na linha i?
linha_teste <- as.numeric(abund_final[i, ])
print("Dados brutos (primeiros 10):")
print(head(linha_teste, 10))

# 2. Quantos NAs?
nas <- sum(is.na(linha_teste))
print(paste("Número de NAs na linha:", nas))

# 3. Quantas espécies presentes (>0)?
presenca_teste <- !is.na(linha_teste) & linha_teste > 0
riqueza_teste <- sum(presenca_teste)
print(paste("Riqueza detectada (S):", riqueza_teste))

# 4. O critério foi atendido?
if(riqueza_teste >= 3) {
  print("SUCESSO: O loop deveria calcular para este sítio!")
  
  # Tenta calcular Pielou
  ab_f_teste <- linha_teste[presenca_teste]
  pielou <- diversity(ab_f_teste) / log(length(ab_f_teste))
  print(paste("Pielou calculado:", round(pielou, 4)))
  
} else {
  print("FALHA: Riqueza menor que 3. O loop vai pular este sítio.")
  print("Possível causa: Interseção de nomes falhou e abund_final tem 0 colunas de espécies?")
  print(paste("Número de colunas em abund_final:", ncol(abund_final)))
}

#================================================================================#
# 3. CÁLCULO DAS MÉTRICAS FUNCIONAIS (ENERGIA vs. FORMA) ------------------------#
#================================================================================#
# 3.1. Limpeza da Biomassa
col_biomassa <- ifelse("dry_weight" %in% colnames(traits_final), "dry_weight", "dry_weigth")
traits_final[[col_biomassa]] <- as.numeric(gsub(",", ".", trimws(as.character(traits_final[[col_biomassa]]))))


# 3.3. Loop de Cálculo TAD/SAD
cols_metricas <- c("SAD_Pielou", "SAD_Slope", 
                   "TAD_Skew_biomass", "TAD_Kurt_biomass",
                   "TAD_Skew_protibia", "TAD_Kurt_protibia",
                   "TAD_Skew_metatibia", "TAD_Kurt_metatibia")
calculos <- data.frame(cod = rownames(abund_final))
calculos[cols_metricas] <- NA


#salvar calculos em .xlsx
write.xlsx(calculos, "calculos_metricas_iniciais.xlsx")

for(i in 1:nrow(abund_final)){
  site_abund <- as.numeric(abund_final[i, ])
  presenca <- site_abund > 0
  ab_f <- site_abund[presenca]
  
  if(sum(presenca) >= 3) {
    # SAD Slope
    ranks <- 1:length(ab_f)
    log_abund <- log(sort(ab_f, decreasing = TRUE))
    modelo_rank <- lm(log_abund ~ ranks)
    calculos$SAD_Slope[i] <- coef(modelo_rank)[2]
    calculos$SAD_Pielou[i] <- diversity(ab_f) / log(length(ab_f))
    
    # TAD Biomassa
    tr_en <- log10(traits_final[[col_biomassa]][presenca])
    pop_en <- rep(tr_en[!is.na(tr_en)], times = round(ab_f[!is.na(tr_en)]))
    if(length(pop_en) >= 3) {
      calculos$TAD_Skew_biomass[i] <- e1071::skewness(pop_en, type=2)
      if(length(pop_en) >= 4) calculos$TAD_Kurt_biomass[i] <- e1071::kurtosis(pop_en, type=2)
    }
    
    # TAD Protibia
    tr_protibia <- traits_final$protibia_area[presenca]
    pop_protibia <- rep(tr_protibia[!is.na(tr_protibia)], times = round(ab_f[!is.na(tr_protibia)]))
    if(length(pop_protibia) >= 3) {
      calculos$TAD_Skew_protibia[i] <- e1071::skewness(pop_protibia, type=2)
      if(length(pop_protibia) >= 4) calculos$TAD_Kurt_protibia[i] <- e1071::kurtosis(pop_protibia, type=2)
    }
    
    # TAD Metatibia
    tr_metatibia <- traits_final$metatibia_length[presenca]
    pop_metatibia <- rep(tr_metatibia[!is.na(tr_metatibia)], times = round(ab_f[!is.na(tr_metatibia)]))
    if(length(pop_metatibia) >= 3) {
      calculos$TAD_Skew_metatibia[i] <- e1071::skewness(pop_metatibia, type=2)
      if(length(pop_metatibia) >= 4) calculos$TAD_Kurt_metatibia[i] <- e1071::kurtosis(pop_metatibia, type=2)
    }
  }
}


# ==============================================================================
# 4. PREPARAÇÃO PARA MODELAGEM
# ==============================================================================

# 1. Garantir que a coluna 'cod' é caractere limpo em ambos
master$cod <- trimws(as.character(master$cod))
calculos$cod <- trimws(as.character(calculos$cod))

# 2. Join Seguro
matriz_analise <- left_join(master, calculos, by = "cod")

# 3. Arrumar Fatores e Escalas
matriz_analise <- matriz_analise %>%
  mutate(uso = factor(uso, levels = c("F", "S", "C", "P"), 
                      labels = c("Forest", "Silviculture", "Coffee", "Pasture"))) %>%
  mutate(
    Floresta_Otimizada = ifelse(uso == "Forest", flo_PLAND_m480, flo_PLAND_m990),
    Focal_Otimizada    = case_when(uso == "Forest" ~ flo_PLAND_m480,
                                   uso == "Pasture" ~ pas_PLAND_m990,
                                   uso == "Coffee" ~ caf_PLAND_m990,
                                   uso == "Silviculture" ~ euc_PLAND_m990),
    Hetero_Otimizada   = ifelse(uso == "Forest", SHDI_m480, SHDI_m990),
    Integ_Otimizada    = ifelse(uso == "Forest", PC1_m480, PC1_m990)
  )

# 4. Calcular PC1 Local (Controle)
vars_locais <- c("temp.max", "umid.min", "canopy_cover", "dens")
matriz_analise$PC1_Local <- NA # Inicializa com NA

if(all(vars_locais %in% colnames(matriz_analise))) {
  completos <- complete.cases(matriz_analise[, vars_locais])
  if(sum(completos) > 3) {
    pca_amb <- prcomp(matriz_analise[completos, vars_locais], scale. = TRUE)
    matriz_analise$PC1_Local[completos] <- pca_amb$x[, 1]
    print("PC1 Local calculado com sucesso.")
  }
} else { warning("Variáveis locais não encontradas!") }






# ==============================================================================
# 5. EXECUÇÃO DOS MODELOS (SELEÇÃO GLOBAL POR MÉTRICA)
# ==============================================================================

vars_paisagem <- c("Floresta_Otimizada", "Focal_Otimizada", "Hetero_Otimizada", "Integ_Otimizada")
lista_aic_global <- list()
lista_vif_global <- list()

print("Iniciando a modelagem global (Todos contra Todos)...")

for(y in cols_metricas) {
  
  # 1. Filtra NAs para TODAS as variáveis envolvidas (Crucial para comparar AICc)
  dados_mod <- matriz_analise %>% 
    filter(!is.na(!!sym(y)), !is.na(PC1_Local)) %>%
    filter(!is.na(Floresta_Otimizada) & !is.na(Focal_Otimizada) & 
             !is.na(Hetero_Otimizada) & !is.na(Integ_Otimizada))
  
  # Se não houver dados suficientes, pula para a próxima métrica
  if(nrow(dados_mod) < 5) next 
  
  # 2. Lista para armazenar TODOS os modelos desta métrica (Y)
  modelos_y <- list()
  
  # Modelo Nulo (m0)
  f_nulo <- as.formula(paste(y, "~ PC1_Local"))
  modelos_y[["m0_Nulo"]] <- glm(f_nulo, data = dados_mod, family = gaussian, na.action = na.fail)
  
  # 3. Loop para criar os modelos aditivos e interativos para cada variável de paisagem
  for(x in vars_paisagem) {
    f_adit <- as.formula(paste(y, "~", x, "+ PC1_Local"))
    f_inte <- as.formula(paste(y, "~", x, "* uso + PC1_Local"))
    
    m1_fit <- glm(f_adit, data = dados_mod, family = gaussian, na.action = na.fail)
    m2_fit <- glm(f_inte, data = dados_mod, family = gaussian, na.action = na.fail)
    
    modelos_y[[paste("m1", x, sep = "_")]] <- m1_fit
    modelos_y[[paste("m2", x, sep = "_")]] <- m2_fit
    
    # Salvar VIF do modelo aditivo
    try({
      vif_res <- car::vif(m1_fit)
      if(is.matrix(vif_res)) val_vif <- vif_res[, 1] else val_vif <- vif_res
      lista_vif_global[[paste(y, x)]] <- data.frame(
        Metrica = y, Variavel_Paisagem = x, Preditor = names(val_vif), VIF = as.numeric(val_vif)
      )
    }, silent = TRUE)
  }
  
  # 4. Seleção de Modelos (Compara os 9 modelos de uma vez)
  tryCatch({
    selecao <- do.call(MuMIn::model.sel, modelos_y)
    
    # Organizar Tabela AIC
    df_aic <- as.data.frame(selecao)
    df_aic$Metrica <- y
    df_aic$Modelo_ID <- rownames(df_aic)
    
    # Extrair qual é a paisagem e o tipo do modelo a partir do nome
    df_aic <- df_aic %>%
      mutate(
        Tipo_Modelo = case_when(grepl("m0", Modelo_ID) ~ "m0", grepl("m1", Modelo_ID) ~ "m1", TRUE ~ "m2"),
        Variavel_Paisagem = case_when(
          grepl("Floresta", Modelo_ID) ~ "Floresta_Otimizada",
          grepl("Focal", Modelo_ID) ~ "Focal_Otimizada",
          grepl("Hetero", Modelo_ID) ~ "Hetero_Otimizada",
          grepl("Integ", Modelo_ID) ~ "Integ_Otimizada",
          TRUE ~ "Nenhuma (Local)"
        )
      ) %>%
      select(Metrica, Modelo_ID, Tipo_Modelo, Variavel_Paisagem, df, logLik, AICc, delta, weight)
    
    lista_aic_global[[y]] <- df_aic
    
  }, error = function(e) message("Erro na seleção para: ", y, " - ", e$message))
}

# 5. Exportar Resultados Globais
tabela_aic_final <- bind_rows(lista_aic_global)
write.xlsx(tabela_aic_final, "Resultados_AICc_Global.xlsx")
write.xlsx(bind_rows(lista_vif_global), "Resultados_VIF_Global.xlsx")

# ==============================================================================
# 6. MANUAL PLOTS IN ENGLISH (MODELS WITH Delta <= 2)
# ==============================================================================

library(ggeffects)
library(ggplot2)
library(cowplot)
library(dplyr)

# Prepare filtered data (Global Filter to ensure consistent N)
data_plot <- matriz_analise %>% 
  filter(!is.na(PC1_Local)) %>%
  filter(!is.na(flo_PLAND_m480) & !is.na(Focal_Otimizada) & 
           !is.na(Hetero_Otimizada) & !is.na(Integ_Otimizada))

# ------------------------------------------------------------------------------
# BLOCK A: MORPHOLOGICAL SKEWNESS (Protibia & Metatibia)
# ------------------------------------------------------------------------------

# A1. Protibia: Integrity (Winner: Weight 0.931)
mod_prot_integ <- glm(TAD_Skew_protibia ~ Integ_Otimizada + PC1_Local, data = data_plot)
pred_prot_integ <- ggpredict(mod_prot_integ, terms = "Integ_Otimizada")

plot_prot_integ <- plot(pred_prot_integ, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "\u0394=0.00 | Weight=0.931",
       x = "Landscape Integrity", y = "TAD Skewness (Protibia)") +
  theme_classic(base_size = 14)

# A2. Metatibia: Focal Habitat (Winner: Weight 0.663)
mod_meta_focal <- glm(TAD_Skew_metatibia ~ Focal_Otimizada + PC1_Local, data = data_plot)
pred_meta_focal <- ggpredict(mod_meta_focal, terms = "Focal_Otimizada")

plot_meta_focal <- plot(pred_meta_focal, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "\u0394=0.00 | Weight=0.663",
       x = "Focal Habitat (%)", y = "TAD Skewness (Metatibia)") +
  theme_classic(base_size = 14)


# ------------------------------------------------------------------------------
# BLOCK B: SAD SLOPE (Abundance Structure - 5 plausible models)
# ------------------------------------------------------------------------------

# B1. Local Environment (m0 - Winner: Weight 0.299)
mod_slope_m0 <- glm(SAD_Slope ~ PC1_Local, data = data_plot)
pred_slope_m0 <- ggpredict(mod_slope_m0, terms = "PC1_Local")

plot_slope_m0 <- plot(pred_slope_m0, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m0 | \u0394=0.00 | Weight=0.299",
       x = "Local Gradient (PC1)", y = "SAD Slope") + theme_classic(base_size = 12)

# B2. Focal Habitat (m1 - Delta 0.46)
mod_slope_focal <- glm(SAD_Slope ~ Focal_Otimizada + PC1_Local, data = data_plot)
pred_slope_focal <- ggpredict(mod_slope_focal, terms = "Focal_Otimizada")

plot_slope_focal <- plot(pred_slope_focal, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=0.46 | Weight=0.238",
       x = "Focal Habitat (%)", y = "SAD Slope") + theme_classic(base_size = 12)

# B3. Heterogeneity (m1 - Delta 1.05)
mod_slope_hetero <- glm(SAD_Slope ~ Hetero_Otimizada + PC1_Local, data = data_plot)
pred_slope_hetero <- ggpredict(mod_slope_hetero, terms = "Hetero_Otimizada")

plot_slope_hetero <- plot(pred_slope_hetero, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=1.05 | Weight=0.177",
       x = "Heterogeneity (SHDI)", y = "SAD Slope") + theme_classic(base_size = 12)

# B4. Forest Cover (m1 - Delta 1.24)
mod_slope_flo <- glm(SAD_Slope ~ Floresta_Otimizada + PC1_Local, data = data_plot)
pred_slope_flo <- ggpredict(mod_slope_flo, terms = "Floresta_Otimizada")

plot_slope_flo <- plot(pred_slope_flo, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=1.24 | Weight=0.161",
       x = "Forest Cover (%)", y = "SAD Slope") + theme_classic(base_size = 12)

# B5. Integrity (m1 - Delta 1.75)
mod_slope_integ <- glm(SAD_Slope ~ Integ_Otimizada + PC1_Local, data = data_plot)
pred_slope_integ <- ggpredict(mod_slope_integ, terms = "Integ_Otimizada")

plot_slope_integ <- plot(pred_slope_integ, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=1.75 | Weight=0.125",
       x = "Landscape Integrity", y = "SAD Slope") + theme_classic(base_size = 12)


# ------------------------------------------------------------------------------
# BLOCK C: KURTOSIS (Biomass)
# ------------------------------------------------------------------------------

# C1. Biomass Kurtosis: Local (m0 - Winner: Weight 0.400)
mod_kurt_bio_m0 <- glm(TAD_Kurt_biomass ~ PC1_Local, data = data_plot)
pred_kurt_bio_m0 <- ggpredict(mod_kurt_bio_m0, terms = "PC1_Local")

plot_kurt_bio_m0 <- plot(pred_kurt_bio_m0, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m0 | \u0394=0.00 | Weight=0.400",
       x = "Local Gradient (PC1)", y = "TAD Kurtosis (Biomass)") + theme_classic()


# C2. Biomass Kurtosis: Focal (m1 - Delta 1.13: weight 0.227) 
mod_kurt_bio_focal <- glm(TAD_Kurt_biomass ~ Focal_Otimizada + PC1_Local, data = data_plot)
pred_kurt_bio_focal <- ggpredict(mod_kurt_bio_focal,terms = "Focal_Otimizada")

plot_kurt_bio_focal <- plot(pred_kurt_bio_focal, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=1.13 | Weight=0.227",
       x = "Focal Habitat (%)", y = "TAD Kurtosis (Biomass)") + theme_classic()


# C3. Biomass Kurtosis: Heterogeneity (m1 - Delta 1.66: weight 0.174)
mod_kurt_bio_hetero <- glm(TAD_Kurt_biomass ~ Hetero_Otimizada + PC1_Local, data = data_plot)
pred_kurt_bio_hetero <- ggpredict(mod_kurt_bio_hetero, terms = "Hetero_Otimizada")

plot_kurt_bio_hetero <- plot(pred_kurt_bio_hetero, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=1.66 | Weight=0.174",
       x = "Heterogeneity (SHDI)", y = "TAD Kurtosis (Biomass)") + theme_classic()


# ------------------------------------------------------------------------------
# BLOCK D: KURTOSIS (Metatibia)
# ------------------------------------------------------------------------------


# D1. Metatibia Kurtosis: Local (m0 - Winner: Weight 0.409)
mod_kurt_meta_m0 <- glm(TAD_Kurt_metatibia ~ PC1_Local, data = data_plot)
pred_kurt_meta_m0 <- ggpredict(mod_kurt_meta_m0, terms = "PC1_Local")

plot_kurt_meta_m0 <- plot(pred_kurt_meta_m0, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m0 | \u0394=0.00 | Weight=0.409",
       x = "Local Gradient (PC1)", y = "TAD Kurtosis (Metatibia)") + theme_classic()


# D2. Metatibia Kurtosis: Focal (m1 - Delta 1.15)
mod_kurt_meta_focal <- glm(TAD_Kurt_metatibia ~ Focal_Otimizada + PC1_Local, data = data_plot)
pred_kurt_meta_focal <- ggpredict(mod_kurt_meta_focal, terms = "Focal_Otimizada")

plot_kurt_meta_focal <- plot(pred_kurt_meta_focal, show_data = TRUE, dot_alpha = 0.5) +
  labs(title = "", subtitle = "m1 | \u0394=1.15 | Weight=0.230",
       x = "Focal Habitat (%)", y = "TAD Kurtosis (Metatibia)") + theme_classic()





# ------------------------------------------------------------------------------
# EXPORTING PANELS (Grids)
# ------------------------------------------------------------------------------

# Panel 1: Morphology Skewness
painel_skew <- plot_grid(plot_prot_integ, plot_meta_focal, ncol = 2, labels = "AUTO")
ggsave("Panel_Morphology_Skewness_EN.png", painel_skew, width = 12, height = 5, dpi = 300)

# Panel 2: SAD Slope (Grid 3x2)
painel_slope <- plot_grid(
  plot_slope_m0, plot_slope_focal, 
  plot_slope_hetero, plot_slope_flo, 
  plot_slope_integ, NULL, 
  ncol = 2, labels = "AUTO", align = "vh"
)
ggsave("Panel_SAD_Slope_EN.png", painel_slope, width = 12, height = 15, dpi = 300)

# Panel 3: Kurtosis (biomass)
painel_kurt <- plot_grid(plot_kurt_bio_m0, plot_kurt_bio_focal, plot_kurt_bio_hetero, ncol = 3, labels = "AUTO")
ggsave("Panel_Kurtosis_EN.png", painel_kurt, width = 12, height = 5, dpi = 300)


# Panel 4: Kurtosis (metatibia)
painel_kurt_meta <- plot_grid(plot_kurt_meta_m0, plot_kurt_meta_focal, ncol = 2, labels = "AUTO")
ggsave("Panel_Kurtosis_Metatibia_EN.png", painel_kurt_meta, width = 12, height = 5, dpi = 300)
