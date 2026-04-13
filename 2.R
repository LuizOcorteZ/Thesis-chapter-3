# ==============================================================================
# SCRIPT DOUTORADO - CAPÍTULO 3 (VERSÃO FINAL CONSOLIDADA - N=35)
# OBJETIVO: Seleção de Modelos (AIC) para Respostas de TAD e SAD
# ABORDAGEM: Hipóteses Concorrentes (Paisagem * Uso + Controle Local via PCA)
# ==============================================================================

# 1. CARREGAR PACOTES ----------------------------------------------------------
pacotes <- c("tidyverse", "vegan", "e1071", "ade4", "openxlsx", 
             "broom", "car", "gridExtra", "grid", "cowplot", "MuMIn") 

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

# 3. CÁLCULO DAS MÉTRICAS (SLOPE SAD + TADS) -----------------------------------

# 3.1. Limpeza da Biomassa
col_biomassa <- ifelse("dry_weight" %in% colnames(traits_final), "dry_weight", "dry_weigth")
traits_final[[col_biomassa]] <- as.numeric(gsub(",", ".", trimws(as.character(traits_final[[col_biomassa]]))))

# 3.2. PCA de Forma (Traits)
traits_shape <- traits_final %>%
  select(where(is.numeric)) %>%
  select(-all_of(col_biomassa), -any_of(c("body_length", "log_dry_weight", "biomass")))

pca_shape <- prcomp(traits_shape, scale. = TRUE, center = TRUE)
traits_final$PC1_Shape <- pca_shape$x[, 1]
traits_final$PC2_Shape <- pca_shape$x[, 2] * -1 

# 3.3. Loop de Cálculo de Métricas por Site
cols_metricas <- c("SAD_Slope", 
                   "TAD_Skew_Biomass", "TAD_Kurt_Biomass",
                   "TAD_Skew_Shape_PC1", "TAD_Kurt_Shape_PC1",
                   "TAD_Skew_Shape_PC2", "TAD_Kurt_Shape_PC2")

calculos <- data.frame(cod = rownames(abund_final))
calculos[cols_metricas] <- NA

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
      calculos$TAD_Skew_Biomass[i] <- e1071::skewness(pop_en, type=2)
      if(length(pop_en) >= 4) calculos$TAD_Kurt_Biomass[i] <- e1071::kurtosis(pop_en, type=2)
    }
    
    # TAD PC1 Shape
    tr_pc1 <- traits_final$PC1_Shape[presenca]
    pop_pc1 <- rep(tr_pc1[!is.na(tr_pc1)], times = round(ab_f[!is.na(tr_pc1)]))
    if(length(pop_pc1) >= 3) {
      calculos$TAD_Skew_Shape_PC1[i] <- e1071::skewness(pop_pc1, type=2)
      if(length(pop_pc1) >= 4) calculos$TAD_Kurt_Shape_PC1[i] <- e1071::kurtosis(pop_pc1, type=2)
    }
    
    # TAD PC2 Shape
    tr_pc2 <- traits_final$PC2_Shape[presenca]
    pop_pc2 <- rep(tr_pc2[!is.na(tr_pc2)], times = round(ab_f[!is.na(tr_pc2)]))
    if(length(pop_pc2) >= 3) {
      calculos$TAD_Skew_Shape_PC2[i] <- e1071::skewness(pop_pc2, type=2)
      if(length(pop_pc2) >= 4) calculos$TAD_Kurt_Shape_PC2[i] <- e1071::kurtosis(pop_pc2, type=2)
    }
  }
}

# 4. PREPARAÇÃO AMBIENTAL E CONTROLE LOCAL -------------------------------------

# 4.1. PCA Ambiental Completo (Sinergia Local) ---------------------------------
# Unindo as variáveis de microclima e estrutura da vegetação
vars_ambientais_locais <- master %>% 
  select(temp.max, umid.med, temp.med, umid.max, temp.min, umid.min, 
         canopy_cover, dens)

# Rodando o PCA único (Essencial usar scale = TRUE)
pca_amb <- prcomp(vars_ambientais_locais, scale. = TRUE, center = TRUE)

# Verificando os pesos para interpretar o eixo
print("--- Pesos (Loadings) do PCA Ambiental Único ---")
print(pca_amb$rotation[, 1]) 

# Salvando o novo controle local (PC1 Ambiental)
master$PC1_Ambiente_Local <- pca_amb$x[, 1]

# DICA: Se os pesos da vegetação forem negativos, multiplique por -1 para 
# facilitar a interpretação (Valores altos = Floresta mais preservada)
master$PC1_Ambiente_Local <- master$PC1_Ambiente_Local * -1

# 4.2. Join e Transformações
matriz_analise <- left_join(master, calculos, by = "cod") %>%
  mutate(uso = factor(uso, levels = c("F", "S", "C", "P"), 
                      labels = c("Forest", "Silviculture", "Coffee", "Pasture")))

# Variáveis de Paisagem Otimizadas
matriz_analise <- matriz_analise %>%
  mutate(
    Floresta_Otimizada = ifelse(uso == "Forest", flo_PLAND_m480, flo_PLAND_m990),
    Focal_Otimizada    = case_when(uso == "Forest" ~ flo_PLAND_m480,
                                   uso == "Pasture" ~ pas_PLAND_m990,
                                   uso == "Coffee" ~ caf_PLAND_m990,
                                   uso == "Silviculture" ~ euc_PLAND_m990),
    Hetero_Otimizada   = ifelse(uso == "Forest", SHDI_m480, SHDI_m990),
    Integ_Otimizada    = ifelse(uso == "Forest", PC1_m480, PC1_m990)
  )

# ==============================================================================
# 5. SELEÇÃO DE MODELOS (AIC) - N=35
# ==============================================================================

hipoteses_paisagem <- list(
  "Habitat_Total"   = "Floresta_Otimizada",
  "Habitat_Focal"   = "Focal_Otimizada",
  "Heterogeneidade" = "Hetero_Otimizada",
  "Integridade"     = "Integ_Otimizada"
)

# CONTROLE LOCAL VIA PCA (Evita inflar erro padrão e colinearidade com microclima)
controles_locais <- "+ PC1_Ambiente_Local"

lista_selecao <- list()
lista_modelos_vencedores <- list()

for(y in cols_metricas) {
  
  # Filtro de dados completos para N=35
  dados_modelagem <- matriz_analise %>% 
    filter(!is.na(!!sym(y)), !is.na(PC1_Ambiente_Local))
  
  modelos_candidatos <- list()
  
  # Modelo Nulo (Referência)
  modelos_candidatos[["Nulo"]] <- glm(as.formula(paste(y, "~ 1")), 
                                      data = dados_modelagem, family = gaussian)
  
  for(nome_h in names(hipoteses_paisagem)) {
    var_x <- hipoteses_paisagem[[nome_h]]
    
    # Modelo com Interação Paisagem*Uso + Controle Local
    f <- as.formula(paste(y, "~", var_x, "* uso", controles_locais))
    
    try({
      mod <- glm(f, data = dados_modelagem, family = gaussian, na.action = na.fail)
      modelos_candidatos[[nome_h]] <- mod
    }, silent = TRUE)
  }
  
  if(length(modelos_candidatos) > 0) {
    selecao <- model.sel(modelos_candidatos)
    df_sel <- as.data.frame(selecao)
    df_sel$Variavel_Resposta <- y
    df_sel$Hipoteses <- rownames(df_sel)
    lista_selecao[[y]] <- df_sel
    
    melhor_nome <- rownames(df_sel)[1]
    lista_modelos_vencedores[[y]] <- modelos_candidatos[[melhor_nome]]
  }
}

# ==============================================================================
# 6. EXTRAÇÃO DOS COEFICIENTES (CORRIGIDO PARA O MODELO NULO)
# ==============================================================================

coefs_vencedores <- list()

# Ordenar a tabela de AIC para garantir que pegamos sempre o vencedor (Delta = 0)
ranking_vencedores <- tabela_aic_final %>% 
  filter(delta == 0) %>% 
  select(Variavel_Resposta, Hipoteses)

for(y in names(lista_modelos_vencedores)) {
  
  mod <- lista_modelos_vencedores[[y]]
  
  # Descobre o nome da hipótese vencedora para essa variável
  nome_hipotese <- ranking_vencedores$Hipoteses[ranking_vencedores$Variavel_Resposta == y]
  
  # --- CONDIÇÃO DE SEGURANÇA ---
  if(nome_hipotese == "Nulo" || length(coef(mod)) == 1) {
    
    # Se o Nulo venceu, cria uma linha informativa sem rodar Anova
    res <- data.frame(
      term = "(Intercepto - Modelo Nulo)",
      sumsq = NA, df = NA, statistic = NA, p.value = NA,
      Variavel_Resposta = y,
      Modelo_Vencedor = "Nulo"
    )
    
  } else {
    
    # Se NÃO for nulo, roda a Anova Tipo III
    tryCatch({
      res <- broom::tidy(car::Anova(mod, type = 3)) %>%
        mutate(
          Variavel_Resposta = y, 
          Modelo_Vencedor = nome_hipotese
        )
    }, error = function(e) {
      # Caso dê erro mesmo não sendo nulo (ex: singularidade)
      res <- data.frame(
        term = "Erro no Ajuste", sumsq = NA, df = NA, statistic = NA, p.value = NA,
        Variavel_Resposta = y, Modelo_Vencedor = nome_hipotese
      )
    })
  }
  
  coefs_vencedores[[y]] <- res
}

# 7. EXPORTAR TABELA FINAL DE ANOVA
tabela_coefs <- bind_rows(coefs_vencedores)

# Remove colunas vazias se houver e organiza
tabela_coefs <- tabela_coefs %>% 
  select(Variavel_Resposta, Modelo_Vencedor, term, statistic, p.value)

write.xlsx(tabela_coefs, "ANOVA_Modelos_Vencedores_Final_Corrigido.xlsx")

print("Tabelas geradas com sucesso! Verifique o arquivo Excel.")
