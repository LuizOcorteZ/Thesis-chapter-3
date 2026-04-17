# ==============================================================================
# SCRIPT DOUTORADO - CAPÍTULO 3
# OBJETIVO: Respostas de TAD e SAD ao Gradiente de Paisagem e Solo
# ==============================================================================

# 1. CONFIGURAÇÕES E PACOTES ---------------------------------------------------
packages <- c("tidyverse", "vegan", "e1071", "ade4", "openxlsx", "MuMIn",
              "broom", "car", "gridExtra", "cowplot", "ggeffects", "ggplot2")

invisible(sapply(packages, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p, quiet = TRUE)
  library(p, character.only = TRUE)
}))

# 2. FUNÇÕES AUXILIARES (HELPER FUNCTIONS) -------------------------------------

#' Calcula métricas de TAD (Trait Abundance Distribution)
calc_tad <- function(trait_vec, ab_vec) {
  # Filtra NAs e alinha abundância
  valid <- !is.na(trait_vec) & !is.na(ab_vec) & ab_vec > 0
  if(sum(ab_vec[valid]) < 3) return(list(skew = NA, kurt = NA))
  
  # Expande a população baseada na abundância (arredondada)
  pop <- rep(trait_vec[valid], times = round(ab_vec[valid]))
  
  res <- list(skew = e1071::skewness(pop, type = 2))
  res$kurt <- if(length(pop) >= 4) e1071::kurtosis(pop, type = 2) else NA
  return(res)
}

#' Função padronizada para plots de modelos
create_model_plot <- function(mod, x_var, x_label, y_label, subtitle = "") {
  pred <- ggpredict(mod, terms = x_var)
  plot(pred, show_data = TRUE, dot_alpha = 0.5) +
    labs(title = "", subtitle = subtitle, x = x_label, y = y_label) +
    theme_classic(base_size = 12)
}

# 3. CARREGAMENTO E LIMPEZA ----------------------------------------------------
master <- read.xlsx("dados_completos_CAP3.xlsx")
abund  <- read.xlsx("abundance_final_itv.xlsx", rowNames = TRUE)
traits <- read.xlsx("traits_final_corrected2.xlsx", rowNames = TRUE)

# Padronização de nomes de colunas/linhas
clean_names <- function(x) gsub("\\.\\.+", ".", gsub("[ _]+", ".", trimws(x)))
colnames(abund) <- clean_names(colnames(abund))
rownames(traits) <- clean_names(rownames(traits))

comuns <- intersect(colnames(abund), rownames(traits))
abund_f  <- abund[, comuns] %>% .[, sort(colnames(.))]
traits_f <- traits[comuns, ] %>% .[order(rownames(.)), ]

stopifnot(all(colnames(abund_f) == rownames(traits_f)))

# Ajuste robusto da Biomassa (corrige o typo 'weigth')
col_bio <- grep("dry_wei", colnames(traits_f), value = TRUE)[1]
traits_f$biomass_log <- log10(as.numeric(gsub(",", ".", traits_f[[col_bio]])))

# 4. CÁLCULO DAS MÉTRICAS ------------------------------------------------------
calculos <- data.frame(cod = rownames(abund_f))

metrics_list <- lapply(1:nrow(abund_f), function(i) {
  ab_i <- as.numeric(abund_f[i, ])
  pres <- ab_i > 0
  
  out <- data.frame(SAD_Pielou = NA, SAD_Slope = NA, 
                    TAD_Skew_bio = NA, TAD_Kurt_bio = NA,
                    TAD_Skew_prot = NA, TAD_Kurt_prot = NA,
                    TAD_Skew_pron = NA, TAD_Kurt_pron = NA,
                    TAD_Skew_met = NA, TAD_Kurt_met = NA)
  
  if(sum(pres) >= 3) {
    # SAD
    ab_f_pres <- ab_i[pres]
    out$SAD_Pielou <- diversity(ab_f_pres) / log(length(ab_f_pres))
    out$SAD_Slope  <- coef(lm(log(sort(ab_f_pres, decreasing = TRUE)) ~ seq_along(ab_f_pres)))[2]
    
    # TADs usando a função auxiliar
    bio  <- calc_tad(traits_f$biomass_log, ab_i)
    prot <- calc_tad(traits_f$protibia_area, ab_i)
    met <- calc_tad(traits_f$metatibia_length, ab_i)
    
    out$TAD_Skew_bio <- bio$skew;  out$TAD_Kurt_bio <- bio$kurt
    out$TAD_Skew_prot <- prot$skew; out$TAD_Kurt_prot <- prot$kurt
    out$TAD_Skew_met <- met$skew; out$TAD_Kurt_met <- met$kurt
  }
  return(out)
})

calculos <- cbind(calculos, bind_rows(metrics_list))

# 5. MODELAGEM E VARIÁVEIS DE SOLO ---------------------------------------------

# ==============================================================================
# 2. TABELA DE REFERÊNCIA DE ESCALAS (SoE Lookup Table)
# ==============================================================================
# ATENÇÃO: Substitua os valores abaixo pelas escalas exatas da sua planilha 
# 'escalas traits.xlsx'. 
# Nota p/ o Café (C): Como Floresta e PC1 não tiveram escala ideal isolada, 
# repetimos a escala do Focal (810) como 'baseline' para ancorar o modelo.

tabela_escalas <- data.frame(
  metrica = c("SAD_Pielou", "SAD_Pielou", "SAD_Pielou",   
              "SAD_Slope","SAD_Slope", "SAD_Slope",
              "TAD_Skew_bio", "TAD_Skew_bio","TAD_Skew_bio",
              "TAD_Kurt_bio", "TAD_Kurt_bio", "TAD_Kurt_bio",
              "TAD_Skew_prot", "TAD_Skew_prot", "TAD_Skew_prot",
              "TAD_Kurt_prot", "TAD_Kurt_prot", "TAD_Kurt_prot",
              "TAD_Skew_pron", "TAD_Skew_pron", "TAD_Skew_pron",
              "TAD_Kurt_pron", "TAD_Kurt_pron", "TAD_Kurt_pron",
              "TAD_Skew_met", "TAD_Skew_met", "TAD_Skew_met",
              "TAD_Kurt_met", "TAD_Kurt_met","TAD_Kurt_met"),
  uso     = c("C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F",
              "C", "S", "F"),
  
  # Escalas em metros para cada variável (Ajuste conforme seus resultados prévios)
  esc_flo  = c(810, 180, 510,
               810, 180, 510,
               810, 390, 510,
               810, 390, 510,
               810, 390, 510,
               810, 390, 510,
               810, 240, 510,
               810, 240, 510,
               810, 240, 510,
               810, 240, 510),
  esc_foc  = c(810, 180, 510,
               810, 180, 510,
               810, 360, 510,
               810, 360, 510,
               810, 360, 510,
               810, 360, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510), # Uso Focal
  esc_pc1  = c(810, 990, 510,
               810, 990, 510,
               810, 570, 510,
               810, 570, 510,
               810, 570, 510,
               810, 570, 510,
               810, 570, 510,
               810, 570, 510,
               810, 570, 510,
               810, 570, 510), # Integridade (810 no C é baseline)
  esc_shdi = c(810, 990, 510,
               810, 990, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510,
               810, 210, 510)  # Heterogeneidade
  )



# Aqui você adiciona as novas colunas de textura (ex: Areia, Argila)
# Certifique-se que elas existam no objeto 'master'
vars_ambientais <- c("Floresta_Otimizada", "Focal_Otimizada", "Hetero_Otimizada", "Integ_Otimizada") # Adicione aqui as novas variáveis!
vars_locais <- c("temp.med", "umid.med", "canopy_cover", "dens")

# Join e Preparação
matriz_analise <- master %>%
  mutate(cod = trimws(as.character(cod))) %>%
  left_join(mutate(calculos, cod = trimws(as.character(cod))), by = "cod") %>%
  mutate(uso = factor(substr(cod, 1, 1))) 


# ==============================================================================
# INSERÇÃO: CÁLCULO DA PCA LOCAL (MICROCLIMA E VEGETAÇÃO)
# ==============================================================================
# Calcula a PCA, padronizando as variáveis (scale.=TRUE) e lidando com NAs
pca_local <- prcomp(matriz_analise[, vars_locais], scale. = TRUE, na.action = na.exclude)

# Anexa o Eixo 1 na matriz principal
matriz_analise$PC_Micro_Veg <- pca_local$x[, 1]
# ==============================================================================

# matriz_analise é o seu data frame contendo todas as variáveis locais e colunas de buffer
matriz_final <- matriz_analise %>%
  filter(uso != "P") %>%                # Exclui a pastagem (N amostral insuficiente)
  mutate(uso = droplevels(factor(uso))) # Remove o nível fantasma para ajustar os Graus de Liberdade



# Exemplo de loop de modelagem AICc simplificado
lista_aic <- list()
metricas_y <- c("SAD_Pielou", "SAD_Slope", "TAD_Skew_bio", "TAD_Kurt_bio",
                "TAD_Skew_prot", "TAD_Kurt_prot", "TAD_Skew_pron", "TAD_Kurt_pron",
                "TAD_Skew_met", "TAD_Kurt_met") # Escolha as métricas para modelar

for(y in metricas_y) {
  
  # A. ISOLAR AS ESCALAS (Extrai a linha do dicionário para a métrica atual)
  ref <- tabela_escalas %>% filter(metrica == y)
  
  # B. CRIAR AS VARIÁVEIS DINÂMICAS DE PAISAGEM E FILTRAR A MATRIZ
  dados_mod <- matriz_final %>%
    mutate(
      # 1. Quantidade de Floresta
      Floresta_SoE = case_when(
        uso == "C" ~ .data[[paste0("flo_PLAND_m", ref$esc_flo[ref$uso=="C"])]],
        uso == "S" ~ .data[[paste0("flo_PLAND_m", ref$esc_flo[ref$uso=="S"])]],
        uso == "F" ~ .data[[paste0("flo_PLAND_m", ref$esc_flo[ref$uso=="F"])]]
      ),
      
      # 2. Quantidade de Uso Focal (Atenção aos prefixos das colunas da sua planilha original)
      Focal_SoE = case_when(
        uso == "C" ~ .data[[paste0("caf_PLAND_m", ref$esc_foc[ref$uso=="C"])]],
        uso == "S" ~ .data[[paste0("euc_PLAND_m", ref$esc_foc[ref$uso=="S"])]],
        uso == "F" ~ .data[[paste0("flo_PLAND_m", ref$esc_foc[ref$uso=="F"])]] 
      ),
      
      # 3. Integridade (PC1 da Paisagem)
      PC1_SoE = case_when(
        uso == "C" ~ .data[[paste0("PC1_m", ref$esc_pc1[ref$uso=="C"])]],
        uso == "S" ~ .data[[paste0("PC1_m", ref$esc_pc1[ref$uso=="S"])]],
        uso == "F" ~ .data[[paste0("PC1_m", ref$esc_pc1[ref$uso=="F"])]]
      ),
      
      # 4. Heterogeneidade (SHDI)
      SHDI_SoE = case_when(
        uso == "C" ~ .data[[paste0("SHDI_m", ref$esc_shdi[ref$uso=="C"])]],
        uso == "S" ~ .data[[paste0("SHDI_m", ref$esc_shdi[ref$uso=="S"])]],
        uso == "F" ~ .data[[paste0("SHDI_m", ref$esc_shdi[ref$uso=="F"])]]
      )
    ) %>%
    
    # A REGRA DE OURO DO AICc: Filtro rigoroso de NAs para manter o 'N' constante
    filter(!is.na(.data[[y]]), 
           !is.na(PC_Micro_Veg),
           !is.na(Floresta_SoE), 
           !is.na(Focal_SoE), 
           !is.na(PC1_SoE), 
           !is.na(SHDI_SoE), 
           !is.na(uso))
  
  # Trava de segurança: Se uma métrica perder muitas áreas por falta de dados, ele avisa
  if(nrow(dados_mod) < 10) {
    warning(paste("Dados insuficientes para", y, "- Pulando para a próxima métrica."))
    next
  }

# ============================================================================
# C. COMPETIÇÃO DE HIPÓTESES (Os 10 Modelos)
# ============================================================================

# Grupo A: Modelos Base
m_null <- glm(as.formula(paste(y, "~ 1")), data = dados_mod)
m_loc  <- glm(as.formula(paste(y, "~ PC_Micro_Veg * uso")), data = dados_mod)

# Grupo B: Hipóteses Puras de Paisagem
m_flo  <- glm(as.formula(paste(y, "~ Floresta_SoE * uso")), data = dados_mod)
m_foc  <- glm(as.formula(paste(y, "~ Focal_SoE * uso")), data = dados_mod)
m_pc1  <- glm(as.formula(paste(y, "~ PC1_SoE * uso")), data = dados_mod)
m_shdi <- glm(as.formula(paste(y, "~ SHDI_SoE * uso")), data = dados_mod)

# Grupo C: Hipóteses Conjuntas (Paisagem + Microclima Local)
m_j_flo  <- glm(as.formula(paste(y, "~ (Floresta_SoE + PC_Micro_Veg) * uso")), data = dados_mod)
m_j_foc  <- glm(as.formula(paste(y, "~ (Focal_SoE + PC_Micro_Veg) * uso")), data = dados_mod)
m_j_pc1  <- glm(as.formula(paste(y, "~ (PC1_SoE + PC_Micro_Veg) * uso")), data = dados_mod)
m_j_shdi <- glm(as.formula(paste(y, "~ (SHDI_SoE + PC_Micro_Veg) * uso")), data = dados_mod)

# ============================================================================
# D. SELEÇÃO E RANQUEAMENTO (Teoria da Informação)
# ============================================================================

# model.sel calcula o AICc, o delta AIC e o peso (weight) de cada hipótese
lista_aic[[y]] <- model.sel(m_null, m_loc, 
                                   m_flo, m_foc, m_pc1, m_shdi,
                                   m_j_flo, m_j_foc, m_j_pc1, m_j_shdi)
}

# ==============================================================================
# 5. INSPEÇÃO DOS RESULTADOS
# ==============================================================================
# Para verificar a tabela de competição para a Assimetria da Biomassa:
print(lista_aic[["TAD_Skew_bio"]])

# Para salvar a tabela do AICc em um arquivo CSV (opcional)
# write.csv(as.data.frame(lista_resultados[["TAD_Skew_bio"]]), "Tabela_AICc_Biomassa.csv")
write.xlsx(bind_rows(lista_aic, .id = "Metrica"), "Tabela_AICc_Competicao_Completa.xlsx")
  
