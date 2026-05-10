# ==============================================================================
# SCRIPT FINAL - CAPÍTULO 3 (Foco no Decaimento Taxonômico e Reshuffling Funcional)
# ==============================================================================

# 1. CONFIGURAÇÕES E PACOTES ---------------------------------------------------
packages <- c("tidyverse", "sads", "MuMIn", "openxlsx", "logistf", "vegan","patchwork","ggplot2")
invisible(sapply(packages, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p, quiet = TRUE)
  library(p, character.only = TRUE)
}))

# Carregamento 
L_raw <- read.xlsx("abundance_final_itv.xlsx", rowNames = TRUE)
L_raw[is.na(L_raw)] <- 0

Q_raw <- read.xlsx("traits_final_corrected2.xlsx", rowNames = TRUE)

R_raw <- read.xlsx("dados_completos_CAP3.xlsx")

# Padronização e Sincronização
padronizar <- function(nome) { gsub("[^a-z0-9]", "", tolower(nome)) }
colnames(L_raw) <- padronizar(colnames(L_raw))
rownames(Q_raw) <- padronizar(rownames(Q_raw))

esp_comuns <- intersect(colnames(L_raw), rownames(Q_raw))
L_f <- L_raw[, esp_comuns]
Q_f <- Q_raw[esp_comuns, ]

# ------------------------------------------------------------------------------
# TRATAMENTO DOS TRAÇOS RESIDUAIS (Garantindo valores positivos para o sads)
# ------------------------------------------------------------------------------
Q_f$protibia_pos <- Q_f$protibia_area - min(Q_f$protibia_area, na.rm = TRUE) + 0.1
Q_f$metatibia_pos <- Q_f$metatibia_length - min(Q_f$metatibia_length, na.rm = TRUE) + 0.1

tracos_alvo <- c(Biomassa = "dry_weigth", 
                 Protibia = "protibia_pos", 
                 Metatibia = "metatibia_pos")

# ------------------------------------------------------------------------------
# FASE 1: DIAGNÓSTICO MULTIDIMENSIONAL (Com critério de Delta AIC >= 2)
# ------------------------------------------------------------------------------
cat("\nIniciando Fase 1: Diagnóstico Multidimensional...\n")
resultados_fase1 <- list()

for (site in rownames(L_f)) {
  abund <- as.numeric(L_f[site, L_f[site,] > 0])
  sp_presentes <- colnames(L_f)[L_f[site,] > 0]
  
  if(length(abund) < 3 || sum(abund) <= 5) next 
  
  abund_int <- round(abund)
  abund_int <- abund_int[abund_int > 0]
  
  # 1. Ajuste da SAD (Taxonomia)
  melhor_sad <- "Erro"
  try({
    m_lnorm <- sads::fitsad(abund_int, "lnorm")
    m_ls    <- sads::fitsad(abund_int, "ls")
    aics    <- c(AIC(m_lnorm), AIC(m_ls))
    delta_aic <- abs(aics[1] - aics[2])
    
    if(delta_aic < 2) {
      melhor_sad <- "Equivalent"
    } else {
      melhor_sad <- ifelse(aics[1] < aics[2], "Log-normal", "Log-series")
    }
  }, silent = TRUE)
  
  site_data <- data.frame(Site = site, Habitat = substr(site, 1, 1), Winner_SAD = melhor_sad)
  
  # 2. Loop para cada Traço
  for (nome_traco in names(tracos_alvo)) {
    coluna_traco <- tracos_alvo[[nome_traco]]
    tracos_vals <- Q_f[sp_presentes, coluna_traco]
    
    valor_tad <- abund * tracos_vals
    tad_int <- round(valor_tad * 1000)
    
    if(length(unique(tad_int)) < 2) {
      tad_int <- tad_int + sample(0:2, length(tad_int), replace = TRUE)
    }
    tad_int <- tad_int[tad_int > 0]
    
    melhor_tad <- "Erro"
    try({
      m_tad_lnorm <- sads::fitsad(tad_int, "lnorm")
      m_tad_ls    <- sads::fitsad(tad_int, "ls")
      aics_tad    <- c(AIC(m_tad_lnorm), AIC(m_tad_ls))
      delta_tad <- abs(aics_tad[1] - aics_tad[2])
      
      if(delta_tad < 2) {
        melhor_tad <- "Equivalent"
      } else {
        melhor_tad <- ifelse(aics_tad[1] < aics_tad[2], "Log-normal", "Log-series")
      }
    }, silent = TRUE)
    
    rs_valor <- suppressWarnings(cor.test(abund, valor_tad, method = "spearman", exact = FALSE)$estimate)
    
    # Salva Y_TAD apenas para tabela descritiva
    site_data[[paste0("Y_TAD_", nome_traco)]] <- ifelse(melhor_tad == "Log-series", 1, 0)
    site_data[[paste0("RS_", nome_traco)]] <- rs_valor
  }
  resultados_fase1[[site]] <- site_data
}

base_diag <- bind_rows(resultados_fase1)
write.xlsx(base_diag, "Diagnostico_Multidimensional_DeltaAIC.xlsx", rowNames = FALSE)

# ------------------------------------------------------------------------------
# FASE 2: PREPARAÇÃO PARA O TORNEIO DE ESCALAS
# ------------------------------------------------------------------------------
base_comum <- base_diag %>% left_join(R_raw, by = c("Site" = "cod"))
escalas <- seq(210, 990, by = 30)

for(esc in escalas) {
  base_comum[[paste0("focal_", esc)]] <- ifelse(base_comum$Habitat == "F", base_comum[[paste0("flo_PLAND_m", esc)]],
                                                ifelse(base_comum$Habitat == "C", base_comum[[paste0("caf_PLAND_m", esc)]],
                                                       ifelse(base_comum$Habitat == "S", base_comum[[paste0("euc_PLAND_m", esc)]],
                                                              base_comum[[paste0("pas_PLAND_m", esc)]])))
}

base_matriz <- base_comum %>% filter(Habitat != "F")
base_floresta <- base_comum %>% filter(Habitat == "F")

# ------------------------------------------------------------------------------
# FASE 3: TORNEIO DE ESCALAS (Scale of Effect) - Foco no Reshuffling (RS)
# ------------------------------------------------------------------------------
cat("\nIniciando Torneio de Escalas...\n")
vencedores_soe <- list()

for (tr in names(tracos_alvo)) {
  df_res_mat_rs <- data.frame()
  df_res_flo_rs <- data.frame()
  
  var_rs <- paste0("RS_", tr)
  
  matriz_limpa_rs <- base_matriz %>% drop_na(all_of(var_rs))
  floresta_limpa_rs <- base_floresta %>% drop_na(all_of(var_rs))
  
  for (esc in escalas) {
    foc <- paste0("focal_", esc)
    flo <- paste0("flo_PLAND_m", esc)
    shd <- paste0("SHDI_m", esc)
    
    # Matriz
    m_foc_rs <- lm(as.formula(paste(var_rs, "~", foc)), data = matriz_limpa_rs)
    m_flo_rs <- lm(as.formula(paste(var_rs, "~", flo)), data = matriz_limpa_rs)
    m_shd_rs <- lm(as.formula(paste(var_rs, "~", shd)), data = matriz_limpa_rs)
    df_res_mat_rs <- rbind(df_res_mat_rs, data.frame(Escala = esc, Focal = AICc(m_foc_rs), Flor = AICc(m_flo_rs), SHDI = AICc(m_shd_rs)))
    
    # Floresta
    m_flo_f <- lm(as.formula(paste(var_rs, "~", flo)), data = floresta_limpa_rs)
    m_shd_f <- lm(as.formula(paste(var_rs, "~", shd)), data = floresta_limpa_rs)
    df_res_flo_rs <- rbind(df_res_flo_rs, data.frame(Escala = esc, Flor = AICc(m_flo_f), SHDI = AICc(m_shd_f)))
  }
  
  vencedores_soe[[tr]] <- list(Matriz_RS = df_res_mat_rs, Floresta_RS = df_res_flo_rs)
}

# Salvar Escalas Vencedoras (Apenas RS)
escalas_vencedoras <- data.frame(
  Traço = names(vencedores_soe),
  Matriz_RS_Focal = sapply(vencedores_soe, function(x) x$Matriz_RS$Escala[which.min(x$Matriz_RS$Focal)]),
  Matriz_RS_Flor  = sapply(vencedores_soe, function(x) x$Matriz_RS$Escala[which.min(x$Matriz_RS$Flor)]),
  Matriz_RS_SHDI  = sapply(vencedores_soe, function(x) x$Matriz_RS$Escala[which.min(x$Matriz_RS$SHDI)]),
  Floresta_RS_Flor = sapply(vencedores_soe, function(x) x$Floresta_RS$Escala[which.min(x$Floresta_RS$Flor)]),
  Floresta_RS_SHDI = sapply(vencedores_soe, function(x) x$Floresta_RS$Escala[which.min(x$Floresta_RS$SHDI)])
)
write.xlsx(escalas_vencedoras, "Escalas_Vencedoras_SOE.xlsx", rowNames = FALSE)


# ------------------------------------------------------------------------------
# FASE 4: MODELAGEM FINAL COM ESCALAS OTIMIZADAS
# ------------------------------------------------------------------------------
vars_micro <- base_comum %>% select(canopy_cover, dens, temp.med, umid.med)
pca_micro <- prcomp(vars_micro, center = TRUE, scale. = TRUE)
base_comum$PC1_micro <- pca_micro$x[, 1]

base_matriz$PC1_micro <- base_comum$PC1_micro[match(base_matriz$Site, base_comum$Site)]
base_floresta$PC1_micro <- base_comum$PC1_micro[match(base_floresta$Site, base_comum$Site)]

get_aicc_firth <- function(model, n) {
  k <- length(model$coefficients)
  ll <- model$loglik[2]
  return(-2*ll + 2*k + (2*k*(k+1))/(n-k-1))
}

n_mat <- nrow(base_matriz)
tabelas_competicao <- list()

# --- A. COMPETIÇÃO DE HIPÓTESES PARA O COLAPSO TAXONÔMICO (Y_SAD) ---
base_matriz$Y_SAD <- ifelse(base_matriz$Winner_SAD == "Log-series", 1, 0)
matriz_sad <- base_matriz %>% drop_na(Y_SAD)

# Escalas fixas (obtidas do seu torneio anterior)
sc_sad <- c(foc = 990, flo = 210, shd = 900)
matriz_sad$foc_opt_sad <- matriz_sad[[paste0("focal_", sc_sad["foc"])]]
matriz_sad$flo_opt_sad <- matriz_sad[[paste0("flo_PLAND_m", sc_sad["flo"])]]
matriz_sad$shd_opt_sad <- matriz_sad[[paste0("SHDI_m", sc_sad["shd"])]]

mSAD_null  <- logistf(Y_SAD ~ 1, data = matriz_sad)
mSAD_flor  <- logistf(Y_SAD ~ flo_opt_sad, data = matriz_sad)
mSAD_flor2 <- logistf(Y_SAD ~ flo_opt_sad + PC1_micro, data = matriz_sad)
mSAD_foc   <- logistf(Y_SAD ~ foc_opt_sad, data = matriz_sad)
mSAD_foc2  <- logistf(Y_SAD ~ foc_opt_sad + PC1_micro, data = matriz_sad)
mSAD_shd   <- logistf(Y_SAD ~ shd_opt_sad, data = matriz_sad)
mSAD_shd2  <- logistf(Y_SAD ~ shd_opt_sad + PC1_micro, data = matriz_sad)
mSAD_loc   <- logistf(Y_SAD ~ PC1_micro, data = matriz_sad)

n_mat_sad <- nrow(matriz_sad)
comp_ysad <- data.frame(
  Hipótese = c("Nulo", "Resgate", "Resgate + Local", "Intensidade", "Intensidade + Local", "Complexidade", "Complexidade + Local", "Local (PC1)"),
  AICc = c(get_aicc_firth(mSAD_null, n_mat_sad), get_aicc_firth(mSAD_flor, n_mat_sad), get_aicc_firth(mSAD_flor2, n_mat_sad), 
           get_aicc_firth(mSAD_foc, n_mat_sad), get_aicc_firth(mSAD_foc2, n_mat_sad), get_aicc_firth(mSAD_shd, n_mat_sad),
           get_aicc_firth(mSAD_shd2, n_mat_sad), get_aicc_firth(mSAD_loc, n_mat_sad))
)
comp_ysad$delta <- comp_ysad$AICc - min(comp_ysad$AICc)
tabelas_competicao[["SAD_Matriz"]] <- comp_ysad %>% arrange(delta)


# --- B. COMPETIÇÃO DE HIPÓTESES PARA RESHUFFLING (RS) ---
escalas_opt <- list(
  Biomassa = list(RS_Matriz = c(foc = 930, flo = 780, shd = 930), RS_Floresta = c(flo = 210, shd = 210)),
  Protibia = list(RS_Matriz = c(foc = 210, flo = 210, shd = 210), RS_Floresta = c(flo = 990, shd = 360)),
  Metatibia = list(RS_Matriz = c(foc = 990, flo = 750, shd = 210), RS_Floresta = c(flo = 210, shd = 210))
)

for (tr in names(escalas_opt)) {
  df_mat <- base_matriz %>% drop_na(all_of(paste0("RS_", tr)))
  df_flo <- base_floresta %>% drop_na(all_of(paste0("RS_", tr)))
  
  var_rs <- paste0("RS_", tr)
  
  # Matriz RS
  sc_M <- escalas_opt[[tr]]$RS_Matriz
  df_mat$foc_opt_rs <- df_mat[[paste0("focal_", sc_M["foc"])]]
  df_mat$flo_opt_rs <- df_mat[[paste0("flo_PLAND_m", sc_M["flo"])]]
  df_mat$shd_opt_rs <- df_mat[[paste0("SHDI_m", sc_M["shd"])]]
  
  mRM_null  <- lm(as.formula(paste(var_rs, "~ 1")), data = df_mat)
  mRM_flor  <- lm(as.formula(paste(var_rs, "~ flo_opt_rs")), data = df_mat)
  mRM_flor2 <- lm(as.formula(paste(var_rs, "~ flo_opt_rs + PC1_micro")), data = df_mat)
  mRM_foc   <- lm(as.formula(paste(var_rs, "~ foc_opt_rs")), data = df_mat)
  mRM_foc2  <- lm(as.formula(paste(var_rs, "~ foc_opt_rs + PC1_micro")), data = df_mat)
  mRM_shd   <- lm(as.formula(paste(var_rs, "~ shd_opt_rs")), data = df_mat)
  mRM_shd2  <- lm(as.formula(paste(var_rs, "~ shd_opt_rs + PC1_micro")), data = df_mat)
  mRM_loc   <- lm(as.formula(paste(var_rs, "~ PC1_micro")), data = df_mat)
  
  comp_rsM <- data.frame(
    Hipótese = c("Nulo", "Resgate", "Resgate + Local", "Intensidade", "Intensidade + Local", "Complexidade", "Complexidade + Local", "Local (PC1)"),
    AICc = c(AICc(mRM_null), AICc(mRM_flor), AICc(mRM_flor2), AICc(mRM_foc), AICc(mRM_foc2), AICc(mRM_shd), AICc(mRM_shd2), AICc(mRM_loc))
  )
  comp_rsM$delta <- comp_rsM$AICc - min(comp_rsM$AICc)
  tabelas_competicao[[paste(tr, "RS_Matriz", sep="_")]] <- comp_rsM %>% arrange(delta)
  
  # Floresta RS
  sc_F <- escalas_opt[[tr]]$RS_Floresta
  df_flo$flo_opt_rs <- df_flo[[paste0("flo_PLAND_m", sc_F["flo"])]]
  df_flo$shd_opt_rs <- df_flo[[paste0("SHDI_m", sc_F["shd"])]]
  
  mRF_null  <- lm(as.formula(paste(var_rs, "~ 1")), data = df_flo)
  mRF_flor  <- lm(as.formula(paste(var_rs, "~ flo_opt_rs")), data = df_flo)
  mRF_flor2 <- lm(as.formula(paste(var_rs, "~ flo_opt_rs + PC1_micro")), data = df_flo)
  mRF_shd   <- lm(as.formula(paste(var_rs, "~ shd_opt_rs")), data = df_flo)
  mRF_shd2  <- lm(as.formula(paste(var_rs, "~ shd_opt_rs + PC1_micro")), data = df_flo)
  mRF_loc   <- lm(as.formula(paste(var_rs, "~ PC1_micro")), data = df_flo)
  
  comp_rsF <- data.frame(
    Hipótese = c("Nulo", "Resgate", "Resgate + Local", "Complexidade", "Complexidade + Local", "Local (PC1)"),
    AICc = c(AICc(mRF_null), AICc(mRF_flor), AICc(mRF_flor2), AICc(mRF_shd), AICc(mRF_shd2), AICc(mRF_loc))
  )
  comp_rsF$delta <- comp_rsF$AICc - min(comp_rsF$AICc)
  tabelas_competicao[[paste(tr, "RS_Floresta", sep="_")]] <- comp_rsF %>% arrange(delta)
}

# ==========================================================================
# EFEITO DA TEXTURA DO SOLO (Somente R_S)
# ==========================================================================
df_solo_completo <- base_comum %>% drop_na(Areia, RS_Biomassa, RS_Protibia, RS_Metatibia)
pca_solo <- prcomp(df_solo_completo %>% select(Areia, Silte, Argila), scale. = TRUE)
df_solo_completo$PC1_solo <- pca_solo$x[, 1]

base_solo_for <- df_solo_completo %>% filter(Habitat == "F")
base_solo_mat <- df_solo_completo %>% filter(Habitat != "F")

tabelas_solo <- list()
for (tr in c("Biomassa", "Protibia", "Metatibia")) {
  var_rs <- paste0("RS_", tr)
  
  # Solo Floresta
  m_null_for <- lm(as.formula(paste(var_rs, "~ 1")), data = base_solo_for)
  m_area_for <- lm(as.formula(paste(var_rs, "~ Areia")), data = base_solo_for)
  comp_for <- data.frame(Variável = var_rs, Hipótese = c("Nulo", "Areia"), AICc = c(AICc(m_null_for), AICc(m_area_for)))
  comp_for$delta <- comp_for$AICc - min(comp_for$AICc)
  tabelas_solo[[paste("Floresta", var_rs, sep="_")]] <- comp_for %>% arrange(delta)
  
  # Solo Matriz
  m_null_mat <- lm(as.formula(paste(var_rs, "~ 1")), data = base_solo_mat)
  m_area_mat <- lm(as.formula(paste(var_rs, "~ Areia")), data = base_solo_mat)
  comp_mat <- data.frame(Variável = var_rs, Hipótese = c("Nulo", "Areia"), AICc = c(AICc(m_null_mat), AICc(m_area_mat)))
  comp_mat$delta <- comp_mat$AICc - min(comp_mat$AICc)
  tabelas_solo[[paste("Matriz", var_rs, sep="_")]] <- comp_mat %>% arrange(delta)
}

# EXPORTAÇÃO DOS MODELOS DE COMPETIÇÃO
df_resultados_paisagem <- bind_rows(tabelas_competicao, .id = "Analise_Variavel")
df_resultados_solo <- bind_rows(tabelas_solo, .id = "Analise_Variavel")
write.xlsx(list("Modelos_Paisagem" = df_resultados_paisagem, "Modelos_Solo" = df_resultados_solo), 
           file = "Resultados_Competicao_Modelos.xlsx", overwrite = TRUE)

# ==========================================================================
# EXTRAÇÃO DA SIGNIFICÂNCIA E INTERVALOS DE CONFIANÇA
# ==========================================================================
vencedores <- list(
  list(ID = "SAD_Matriz",        Base = "base_matriz",  Y = "Y_SAD",        X = "focal_990",  Tipo = "Firth"),
  list(ID = "RS_Bio_Matriz",     Base = "base_matriz",  Y = "RS_Biomassa",  X = "focal_930",  Tipo = "Linear"),
  list(ID = "RS_Bio_Floresta",   Base = "base_floresta",Y = "RS_Biomassa",  X = "SHDI_m210",  Tipo = "Linear"),
  list(ID = "RS_Prot_Matriz",    Base = "base_matriz",  Y = "RS_Protibia",  X = "SHDI_m210",  Tipo = "Linear"),
  list(ID = "RS_Meta_Floresta",  Base = "base_floresta",Y = "RS_Metatibia", X = "SHDI_m210",  Tipo = "Linear"),
  list(ID = "Solo_RS_Prot_Flo",  Base = "base_solo_for",Y = "RS_Protibia",  X = "Areia",      Tipo = "Linear")
)

tabela_final <- data.frame()
for (m in vencedores) {
  dados <- get(m$Base) %>% drop_na(all_of(c(m$Y, m$X)))
  formula_mod <- as.formula(paste(m$Y, "~", m$X))
  
  if (m$Tipo == "Linear") {
    mod <- lm(formula_mod, data = dados)
    sum_mod <- summary(mod)
    ic <- confint(mod)
    resumo <- data.frame(
      Analise = m$ID, Variavel_X = m$X, Estimate = sum_mod$coefficients[2, 1], Std_Error = sum_mod$coefficients[2, 2],
      P_Value = sum_mod$coefficients[2, 4], IC_2.5 = ic[2, 1], IC_97.5 = ic[2, 2], R2_Adj = sum_mod$adj.r.squared
    )
  } else {
    mod <- logistf(formula_mod, data = dados)
    sum_mod <- summary(mod)
    resumo <- data.frame(
      Analise = m$ID, Variavel_X = m$X, Estimate = mod$coefficients[2], Std_Error = sqrt(diag(mod$var))[2],
      P_Value = sum_mod$prob[2], IC_2.5 = mod$ci.lower[2], IC_97.5 = mod$ci.upper[2], R2_Adj = NA
    )
  }
  tabela_final <- rbind(tabela_final, resumo)
}

tabela_final$Significativo <- ifelse(tabela_final$IC_2.5 * tabela_final$IC_97.5 > 0, "SIM", "NÃO")
write.xlsx(tabela_final, "Significancia_e_Intervalos_Confianca_Limpo.xlsx", overwrite = TRUE)
cat("\nScript Final concluído! Tabelas prontas.\n")
