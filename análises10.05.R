# ==============================================================================
# SCRIPT FINAL - CAPÍTULO 3 (Foco no Decaimento Taxonômico e Reshuffling Funcional)
# ==============================================================================

# 1. CONFIGURAÇÕES E PACOTES ---------------------------------------------------
# Definindo semente para garantir reprodutibilidade (crucial para o jitter das TADs)
set.seed(1234)

packages <- c("tidyverse", "sads", "MuMIn", "openxlsx", "logistf", "vegan", "patchwork", "ggplot2", "here")
invisible(sapply(packages, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p, quiet = TRUE)
  library(p, character.only = TRUE)
}))

# ------------------------------------------------------------------------------
# IMPORTAÇÃO E SINCRONIZAÇÃO
# ------------------------------------------------------------------------------
message("Carregando e sincronizando dados...")

# Uso do here() para caminhos relativos robustos (evita erros em computadores diferentes)
L_raw <- read.xlsx(here::here("abundance_final_itv.xlsx"), rowNames = TRUE)
L_raw[is.na(L_raw)] <- 0

Q_raw <- read.xlsx(here::here("traits_final_corrected2.xlsx"), rowNames = TRUE)
R_raw <- read.xlsx(here::here("dados_completos_CAP3.xlsx"))

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
message("Iniciando Fase 1: Diagnóstico Multidimensional...")

# Função auxiliar para ajustar as distribuições e escolher a vencedora
ajustar_sad <- function(abundances) {
  tryCatch({
    m_lnorm <- sads::fitsad(abundances, "lnorm")
    m_ls    <- sads::fitsad(abundances, "ls")
    aics    <- c(AIC(m_lnorm), AIC(m_ls))
    delta_aic <- abs(aics[1] - aics[2])
    
    if (delta_aic < 2) return("Equivalent")
    return(ifelse(aics[1] < aics[2], "Log-normal", "Log-series"))
  }, error = function(e) return("Erro"))
}

# Função principal para processar os dados de um único local (Site)
processar_site <- function(site, L, Q, tracos) {
  abund <- as.numeric(L[site, L[site,] > 0])
  sp_presentes <- colnames(L)[L[site,] > 0]
  
  # Regra de corte para exclusão de sites pobres
  if(length(abund) < 3 || sum(abund) <= 5) return(NULL) 
  
  abund_int <- round(abund)
  abund_int <- abund_int[abund_int > 0]
  
  # 1. Ajuste da SAD (Taxonomia)
  melhor_sad <- ajustar_sad(abund_int)
  site_data <- data.frame(Site = site, Habitat = substr(site, 1, 1), Winner_SAD = melhor_sad)
  
  # 2. Ajuste para cada Traço Funcional (TAD)
  for (nome_traco in names(tracos)) {
    coluna_traco <- tracos[[nome_traco]]
    tracos_vals <- Q[sp_presentes, coluna_traco]
    
    valor_tad <- abund * tracos_vals
    tad_int <- round(valor_tad * 1000)
    
    # Adição de ruído (jitter) para desempate
    if(length(unique(tad_int)) < 2) {
      tad_int <- tad_int + sample(0:2, length(tad_int), replace = TRUE)
    }
    tad_int <- tad_int[tad_int > 0]
    
    melhor_tad <- ajustar_sad(tad_int)
    rs_valor <- suppressWarnings(cor.test(abund, valor_tad, method = "spearman", exact = FALSE)$estimate)
    
    # Salva Y_TAD para tabela descritiva e Reshuffling Funcional (RS)
    site_data[[paste0("Y_TAD_", nome_traco)]] <- ifelse(melhor_tad == "Log-series", 1, 0)
    site_data[[paste0("RS_", nome_traco)]] <- rs_valor
  }
  
  return(site_data)
}

# Aplicando a função a todos os sites de forma vetorizada
lista_resultados <- lapply(rownames(L_f), processar_site, L = L_f, Q = Q_f, tracos = tracos_alvo)

# Consolidando a lista em um único dataframe de forma eficiente
base_diag <- bind_rows(lista_resultados)

write.xlsx(base_diag, here::here("Diagnostico_Multidimensional_DeltaAIC.xlsx"), rowNames = FALSE)
message(sprintf("Fase 1 concluída! %d sites processados com sucesso e exportados.", nrow(base_diag)))



# ------------------------------------------------------------------------------
# FASE 2: PREPARAÇÃO PARA O TORNEIO DE ESCALAS
# ------------------------------------------------------------------------------
message("\nIniciando Fase 2: Preparação das Escalas...")

base_comum <- base_diag %>% left_join(R_raw, by = c("Site" = "cod"))
escalas <- seq(210, 990, by = 30)

# Uso de case_when para tornar as regras de Habitat legíveis e seguras
for(esc in escalas) {
  col_flo <- paste0("flo_PLAND_m", esc)
  col_caf <- paste0("caf_PLAND_m", esc)
  col_euc <- paste0("euc_PLAND_m", esc)
  col_pas <- paste0("pas_PLAND_m", esc)
  col_foc <- paste0("focal_", esc)
  
  base_comum <- base_comum %>%
    mutate(!!col_foc := case_when(
      Habitat == "F" ~ .data[[col_flo]],
      Habitat == "C" ~ .data[[col_caf]],
      Habitat == "S" ~ .data[[col_euc]],
      TRUE           ~ .data[[col_pas]]
    ))
}

base_matriz <- base_comum %>% filter(Habitat != "F")
base_floresta <- base_comum %>% filter(Habitat == "F")

# ------------------------------------------------------------------------------
# FASE 3: TORNEIO DE ESCALAS (Scale of Effect) - Foco no Reshuffling (RS)
# ------------------------------------------------------------------------------
message("\nIniciando Fase 3: Torneio de Escalas (Scale of Effect)...")
vencedores_soe <- list()

for (tr in names(tracos_alvo)) {
  var_rs <- paste0("RS_", tr)
  
  matriz_limpa_rs <- base_matriz %>% drop_na(all_of(var_rs))
  floresta_limpa_rs <- base_floresta %>% drop_na(all_of(var_rs))
  
  # Informa o tamanho amostral após remover os NAs (excelente para a Tese)
  message(sprintf("Calculando SoE para %s | N_Matriz = %d, N_Floresta = %d", 
                  tr, nrow(matriz_limpa_rs), nrow(floresta_limpa_rs)))
  
  # Usando lapply em vez de rbind no loop: ganho gigantesco de performance
  resultados_escala <- lapply(escalas, function(esc) {
    foc <- paste0("focal_", esc)
    flo <- paste0("flo_PLAND_m", esc)
    shd <- paste0("SHDI_m", esc)
    
    # Modelos Matriz
    m_foc_rs <- lm(as.formula(paste(var_rs, "~", foc)), data = matriz_limpa_rs)
    m_flo_rs <- lm(as.formula(paste(var_rs, "~", flo)), data = matriz_limpa_rs)
    m_shd_rs <- lm(as.formula(paste(var_rs, "~", shd)), data = matriz_limpa_rs)
    
    df_mat <- data.frame(Escala = esc, 
                         Focal = AICc(m_foc_rs), 
                         Flor = AICc(m_flo_rs), 
                         SHDI = AICc(m_shd_rs))
    
    # Modelos Floresta
    m_flo_f <- lm(as.formula(paste(var_rs, "~", flo)), data = floresta_limpa_rs)
    m_shd_f <- lm(as.formula(paste(var_rs, "~", shd)), data = floresta_limpa_rs)
    
    df_flo <- data.frame(Escala = esc, 
                         Flor = AICc(m_flo_f), 
                         SHDI = AICc(m_shd_f))
    
    return(list(mat = df_mat, flo = df_flo))
  })
  
  # Extraindo e unindo as tabelas de resultados de forma otimizada
  df_res_mat_rs <- bind_rows(lapply(resultados_escala, `[[`, "mat"))
  df_res_flo_rs <- bind_rows(lapply(resultados_escala, `[[`, "flo"))
  
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

write.xlsx(escalas_vencedoras, here::here("Escalas_Vencedoras_SOE.xlsx"), rowNames = FALSE)
message("Torneio de Escalas concluído! Arquivo 'Escalas_Vencedoras_SOE.xlsx' gerado.")



# ------------------------------------------------------------------------------
# FASE 4: MODELAGEM FINAL E IMPORTÂNCIA RELATIVA (SEM FDR)
# ------------------------------------------------------------------------------
message("\nIniciando Fase 4: Modelagem Final e Importância Relativa...")

# 4.1 PREPARAÇÃO DOS DADOS (Microclima)
vars_micro <- base_comum %>% select(canopy_cover, dens, temp.med, umid.med) %>% drop_na()
pca_micro <- prcomp(vars_micro, center = TRUE, scale. = TRUE)
base_comum$PC1_micro <- pca_micro$x[match(rownames(base_comum), rownames(vars_micro)), 1]

base_matriz$PC1_micro <- base_comum$PC1_micro[match(base_matriz$Site, base_comum$Site)]
base_floresta$PC1_micro <- base_comum$PC1_micro[match(base_floresta$Site, base_comum$Site)]

get_aicc_firth <- function(model, n) {
  k <- length(model$coefficients)
  ll <- model$loglik[2]
  return(-2*ll + 2*k + (2*k*(k+1))/(n-k-1))
}

# 4.2 CÁLCULO DE IMPORTÂNCIA RELATIVA
calcular_importancia <- function(tabela_comp) {
  melhores <- tabela_comp %>% filter(delta <= 2)
  if (nrow(melhores) <= 1 || "Nulo" %in% melhores$Hipótese) return(NA)
  
  melhores$peso_akaike <- exp(-0.5 * melhores$delta) / sum(exp(-0.5 * melhores$delta))
  vars_map <- list(
    "Resgate" = "flor", "Resgate + Local" = c("flor", "PC1_micro"),
    "Intensidade" = "foc", "Intensidade + Local" = c("foc", "PC1_micro"),
    "Complexidade" = "shd", "Complexidade + Local" = c("shd", "PC1_micro"),
    "Local (PC1)" = "PC1_micro", "Areia" = "Areia"
  )
  
  importancias <- list()
  for (hip in melhores$Hipótese) {
    peso <- melhores$peso_akaike[melhores$Hipótese == hip]
    for (v in vars_map[[hip]]) importancias[[v]] <- sum(importancias[[v]] %||% 0, peso)
  }
  return(paste(names(importancias), round(unlist(importancias), 3), sep = ": ", collapse = " | "))
}

tabelas_competicao <- list()
importancia_resultados <- data.frame(Analise = character(), Variaveis_Importantes = character())
resultados_finais <- list()

# ------------------------------------------------------------------------------
# A. COMPETIÇÃO E EXTRAÇÃO - COLAPSO TAXONÔMICO (Y_SAD)
# ------------------------------------------------------------------------------
base_matriz$Y_SAD <- ifelse(base_matriz$Winner_SAD == "Log-series", 1, 0)
matriz_sad <- base_matriz %>% drop_na(Y_SAD, PC1_micro)

sc_sad <- c(foc = 990, flo = 210, shd = 900) 
matriz_sad$foc_opt <- matriz_sad[[paste0("focal_", sc_sad["foc"])]]
matriz_sad$flo_opt <- matriz_sad[[paste0("flo_PLAND_m", sc_sad["flo"])]]
matriz_sad$shd_opt <- matriz_sad[[paste0("SHDI_m", sc_sad["shd"])]]

mSAD_null  <- logistf(Y_SAD ~ 1, data = matriz_sad)
mSAD_flor  <- logistf(Y_SAD ~ flo_opt, data = matriz_sad)
mSAD_flor2 <- logistf(Y_SAD ~ flo_opt + PC1_micro, data = matriz_sad)
mSAD_foc   <- logistf(Y_SAD ~ foc_opt, data = matriz_sad)
mSAD_foc2  <- logistf(Y_SAD ~ foc_opt + PC1_micro, data = matriz_sad)
mSAD_shd   <- logistf(Y_SAD ~ shd_opt, data = matriz_sad)
mSAD_shd2  <- logistf(Y_SAD ~ shd_opt + PC1_micro, data = matriz_sad)
mSAD_loc   <- logistf(Y_SAD ~ PC1_micro, data = matriz_sad)

n_sad <- nrow(matriz_sad)
comp_ysad <- data.frame(
  Hipótese = c("Nulo", "Resgate", "Resgate + Local", "Intensidade", "Intensidade + Local", "Complexidade", "Complexidade + Local", "Local (PC1)"),
  AICc = c(get_aicc_firth(mSAD_null, n_sad), get_aicc_firth(mSAD_flor, n_sad), get_aicc_firth(mSAD_flor2, n_sad), 
           get_aicc_firth(mSAD_foc, n_sad), get_aicc_firth(mSAD_foc2, n_sad), get_aicc_firth(mSAD_shd, n_sad),
           get_aicc_firth(mSAD_shd2, n_sad), get_aicc_firth(mSAD_loc, n_sad))
) %>% mutate(delta = AICc - min(AICc)) %>% arrange(delta)

tabelas_competicao[["SAD_Matriz"]] <- comp_ysad
importancia_resultados <- rbind(importancia_resultados, data.frame(Analise = "SAD_Matriz", Variaveis_Importantes = calcular_importancia(comp_ysad)))

# Extração SAD (Exemplo usando Focal 990)
sum_sad <- summary(mSAD_foc)
resultados_finais[["SAD_Matriz"]] <- data.frame(
  Analise = "SAD_Matriz", Variavel_X = "focal_990", 
  Estimate = mSAD_foc$coefficients[2], Std_Error = sqrt(diag(mSAD_foc$var))[2],
  P_Value = sum_sad$prob[2], IC_2.5 = mSAD_foc$ci.lower[2], IC_97.5 = mSAD_foc$ci.upper[2]
)

# ==============================================================================
# B. COMPETIÇÃO E EXTRAÇÃO (LRT) - RESHUFFLING FUNCIONAL (RS)
# ==============================================================================
escalas_opt <- list(
  Biomassa  = list(RS_Matriz = c(foc = 930, flo = 780, shd = 930), RS_Floresta = c(flo = 210, shd = 210)),
  Protibia  = list(RS_Matriz = c(foc = 210, flo = 210, shd = 210), RS_Floresta = c(flo = 990, shd = 360)),
  Metatibia = list(RS_Matriz = c(foc = 990, flo = 750, shd = 210), RS_Floresta = c(flo = 210, shd = 210))
)

for (tr in names(escalas_opt)) {
  var_rs <- paste0("RS_", tr)
  df_mat <- base_matriz %>% drop_na(all_of(var_rs), PC1_micro)
  df_flo <- base_floresta %>% drop_na(all_of(var_rs), PC1_micro)
  
  # --- 1. MATRIZ (RS) ---
  sc_M <- escalas_opt[[tr]]$RS_Matriz
  df_mat$foc_opt <- df_mat[[paste0("focal_", sc_M["foc"])]]
  df_mat$flo_opt <- df_mat[[paste0("flo_PLAND_m", sc_M["flo"])]]
  df_mat$shd_opt <- df_mat[[paste0("SHDI_m", sc_M["shd"])]]
  
  mRM_null  <- lm(as.formula(paste(var_rs, "~ 1")), data = df_mat)
  mRM_flor  <- lm(as.formula(paste(var_rs, "~ flo_opt")), data = df_mat)
  mRM_flor2 <- lm(as.formula(paste(var_rs, "~ flo_opt + PC1_micro")), data = df_mat)
  mRM_foc   <- lm(as.formula(paste(var_rs, "~ foc_opt")), data = df_mat)
  mRM_foc2  <- lm(as.formula(paste(var_rs, "~ foc_opt + PC1_micro")), data = df_mat)
  mRM_shd   <- lm(as.formula(paste(var_rs, "~ shd_opt")), data = df_mat)
  mRM_shd2  <- lm(as.formula(paste(var_rs, "~ shd_opt + PC1_micro")), data = df_mat)
  mRM_loc   <- lm(as.formula(paste(var_rs, "~ PC1_micro")), data = df_mat)
  
  comp_rsM <- data.frame(
    Hipótese = c("Nulo", "Resgate", "Resgate + Local", "Intensidade", "Intensidade + Local", "Complexidade", "Complexidade + Local", "Local (PC1)"),
    AICc = c(AICc(mRM_null), AICc(mRM_flor), AICc(mRM_flor2), AICc(mRM_foc), AICc(mRM_foc2), AICc(mRM_shd), AICc(mRM_shd2), AICc(mRM_loc))
  ) %>% mutate(delta = AICc - min(AICc)) %>% arrange(delta)
  
  nome_mat <- paste(tr, "RS_Matriz", sep="_")
  tabelas_competicao[[nome_mat]] <- comp_rsM
  importancia_resultados <- rbind(importancia_resultados, data.frame(Analise = nome_mat, Variaveis_Importantes = calcular_importancia(comp_rsM)))
  
  # Extração LRT: Comparando Modelo Vencedor (ex: Focal) com o Nulo
  lrt_mat <- anova(mRM_null, mRM_foc, test = "Chisq") 
  p_lrt_mat <- lrt_mat$`Pr(>Chi)`[2]
  
  melhor_mod_mat <- summary(mRM_foc) 
  ic_mat <- confint(mRM_foc)
  
  resultados_finais[[nome_mat]] <- data.frame(
    Analise = nome_mat, Variavel_X = paste0("focal_", sc_M["foc"]),
    Estimate = melhor_mod_mat$coefficients[2, 1], Std_Error = melhor_mod_mat$coefficients[2, 2],
    P_Value_LRT = p_lrt_mat, IC_2.5 = ic_mat[2, 1], IC_97.5 = ic_mat[2, 2]
  )
  
  # --- 2. FLORESTA (RS) ---
  sc_F <- escalas_opt[[tr]]$RS_Floresta
  df_flo$flo_opt <- df_flo[[paste0("flo_PLAND_m", sc_F["flo"])]]
  df_flo$shd_opt <- df_flo[[paste0("SHDI_m", sc_F["shd"])]]
  
  mRF_null  <- lm(as.formula(paste(var_rs, "~ 1")), data = df_flo)
  mRF_flor  <- lm(as.formula(paste(var_rs, "~ flo_opt")), data = df_flo)
  mRF_flor2 <- lm(as.formula(paste(var_rs, "~ flo_opt + PC1_micro")), data = df_flo)
  mRF_shd   <- lm(as.formula(paste(var_rs, "~ shd_opt")), data = df_flo)
  mRF_shd2  <- lm(as.formula(paste(var_rs, "~ shd_opt + PC1_micro")), data = df_flo)
  mRF_loc   <- lm(as.formula(paste(var_rs, "~ PC1_micro")), data = df_flo)
  
  comp_rsF <- data.frame(
    Hipótese = c("Nulo", "Resgate", "Resgate + Local", "Complexidade", "Complexidade + Local", "Local (PC1)"),
    AICc = c(AICc(mRF_null), AICc(mRF_flor), AICc(mRF_flor2), AICc(mRF_shd), AICc(mRF_shd2), AICc(mRF_loc))
  ) %>% mutate(delta = AICc - min(AICc)) %>% arrange(delta)
  
  nome_flo <- paste(tr, "RS_Floresta", sep="_")
  tabelas_competicao[[nome_flo]] <- comp_rsF
  importancia_resultados <- rbind(importancia_resultados, data.frame(Analise = nome_flo, Variaveis_Importantes = calcular_importancia(comp_rsF)))
  
  # Extração LRT: Comparando Modelo Vencedor (ex: SHDI) com o Nulo
  lrt_flo <- anova(mRF_null, mRF_shd, test = "Chisq")
  p_lrt_flo <- lrt_flo$`Pr(>Chi)`[2]
  
  melhor_mod_flo <- summary(mRF_shd) 
  ic_flo <- confint(mRF_shd)
  
  resultados_finais[[nome_flo]] <- data.frame(
    Analise = nome_flo, Variavel_X = paste0("SHDI_m", sc_F["shd"]),
    Estimate = melhor_mod_flo$coefficients[2, 1], Std_Error = melhor_mod_flo$coefficients[2, 2],
    P_Value_LRT = p_lrt_flo, IC_2.5 = ic_flo[2, 1], IC_97.5 = ic_flo[2, 2]
  )
}

# ==============================================================================
# C. COMPETIÇÃO E EXTRAÇÃO (LRT) - EFEITO DA TEXTURA DO SOLO (RS)
# ==============================================================================

# 1. DEFINIR A BASE DE SOLO (Unindo os resultados de RS com os dados originais)
df_solo_completo <- base_comum %>% drop_na(RS_Biomassa, RS_Protibia, RS_Metatibia)

# 2. SEPARAR POR HABITAT (E remover os NAs da coluna 'Areia' apenas aqui, de forma segura)
base_solo_for <- df_solo_completo %>% filter(Habitat == "F")
base_solo_mat <- df_solo_completo %>% filter(Habitat != "F")

for (tr in c("Biomassa", "Protibia", "Metatibia")) {
  var_rs <- paste0("RS_", tr)
  
  # --- 1. MATRIZ (SOLO) ---
  # Cria um dataframe temporário sem nenhum NA nas variáveis específicas do modelo atual
  df_mat_limpo <- base_solo_mat %>% drop_na(all_of(c(var_rs, "Areia")))
  
  if(nrow(df_mat_limpo) > 0) {
    mSM_null  <- lm(as.formula(paste(var_rs, "~ 1")), data = df_mat_limpo)
    mSM_areia <- lm(as.formula(paste(var_rs, "~ Areia")), data = df_mat_limpo)
    
    # Teste LRT Solo Matriz
    lrt_soloM <- anova(mSM_null, mSM_areia, test = "Chisq")
    p_lrt_soloM <- lrt_soloM$`Pr(>Chi)`[2]
    
    sum_soloM <- summary(mSM_areia)
    ic_soloM <- confint(mSM_areia)
    
    nome_soloM <- paste("Solo_Matriz", tr, sep="_")
    resultados_finais[[nome_soloM]] <- data.frame(
      Analise = nome_soloM, Variavel_X = "Areia",
      Estimate = sum_soloM$coefficients[2, 1], Std_Error = sum_soloM$coefficients[2, 2],
      P_Value_LRT = p_lrt_soloM, IC_2.5 = ic_soloM[2, 1], IC_97.5 = ic_soloM[2, 2]
    )
  }
  
  # --- 2. FLORESTA (SOLO) ---
  # Cria um dataframe temporário sem nenhum NA nas variáveis específicas do modelo atual
  df_for_limpo <- base_solo_for %>% drop_na(all_of(c(var_rs, "Areia")))
  
  if(nrow(df_for_limpo) > 0) {
    mSF_null  <- lm(as.formula(paste(var_rs, "~ 1")), data = df_for_limpo)
    mSF_areia <- lm(as.formula(paste(var_rs, "~ Areia")), data = df_for_limpo)
    
    # Teste LRT Solo Floresta
    lrt_soloF <- anova(mSF_null, mSF_areia, test = "Chisq")
    p_lrt_soloF <- lrt_soloF$`Pr(>Chi)`[2]
    
    sum_soloF <- summary(mSF_areia)
    ic_soloF <- confint(mSF_areia)
    
    nome_soloF <- paste("Solo_Floresta", tr, sep="_")
    resultados_finais[[nome_soloF]] <- data.frame(
      Analise = nome_soloF, Variavel_X = "Areia",
      Estimate = sum_soloF$coefficients[2, 1], Std_Error = sum_soloF$coefficients[2, 2],
      P_Value_LRT = p_lrt_soloF, IC_2.5 = ic_soloF[2, 1], IC_97.5 = ic_soloF[2, 2]
    )
  }
}

# ==============================================================================
# CONSOLIDAÇÃO FINAL COM P-VALOR DO LRT
# ==============================================================================
tabela_completa <- bind_rows(resultados_finais)

# Adiciona a SAD (que já tinha sido rodada no bloco anterior e cujo P_Value já é LRT/PLRT)
# Apenas renomeamos a coluna para P_Value_LRT para parear com as outras
resultados_finais[["SAD_Matriz"]] <- resultados_finais[["SAD_Matriz"]] %>% 
  rename(P_Value_LRT = P_Value)

tabela_completa <- bind_rows(resultados_finais)

# Critério de Significância Atualizado: IC de 95% + LRT p < 0.05
tabela_completa$Significativo <- ifelse(
  tabela_completa$IC_2.5 * tabela_completa$IC_97.5 > 0 & tabela_completa$P_Value_LRT < 0.05, 
  "SIM", "NÃO"
)

tabela_completa <- tabela_completa %>% 
  left_join(importancia_resultados, by = "Analise") %>%
  distinct(Analise, .keep_all = TRUE)

# Exportar
lista_exportacao <- list(
  "Resumo_Modelos" = tabela_completa,
  "Competicao_Paisagem" = bind_rows(tabelas_competicao, .id = "Analise")
)
write.xlsx(lista_exportacao, here::here("Resultados_Finais_LRT.xlsx"), overwrite = TRUE)
message("Análise finalizada! Arquivo 'Resultados_Finais_LRT.xlsx' gerado com os resultados consolidados.")





# ==============================================================================
# GERAÇÃO DOS GRÁFICOS PARA MODELOS SIGNIFICATIVOS
# ==============================================================================
message("\nGerando visualizações (scatterplots) para os modelos significativos...")

# 1. Biomassa na Matriz (~ focal_930)
p1 <- ggplot(base_matriz %>% drop_na(RS_Biomassa, focal_930), aes(x = focal_930, y = RS_Biomassa)) +
  geom_point(alpha = 0.7, size = 3, color = "#E69F00") +
  geom_smooth(method = "lm", color = "black", fill = "gray80") +
  theme_classic() +
  labs(title = "A) Matrix Biomass",
       x = " Landscape Intensification (Focal 930m)",
       y = "Functional Reshuffling (Biomass)")

# 2. Biomassa na Floresta (~ SHDI_m210)
p2 <- ggplot(base_floresta %>% drop_na(RS_Biomassa, SHDI_m210), aes(x = SHDI_m210, y = RS_Biomassa)) +
  geom_point(alpha = 0.7, size = 3, color = "#56B4E9") +
  geom_smooth(method = "lm", color = "black", fill = "gray80") +
  theme_classic() +
  labs(title = "B) Forest biomass",
       x = "Landscape Heterogeneity (SHDI 210m)",
       y = "Functional Reshuffling (Biomass)")

# 3. Metatíbia na Floresta (~ SHDI_m210)
p3 <- ggplot(base_floresta %>% drop_na(RS_Metatibia, SHDI_m210), aes(x = SHDI_m210, y = RS_Metatibia)) +
  geom_point(alpha = 0.7, size = 3, color = "#009E73") +
  geom_smooth(method = "lm", color = "black", fill = "gray80") +
  theme_classic() +
  labs(title = "C) Forest Metatibia",
       x = "Landscape Heterogeneity (SHDI 210m)",
       y = "Functional Reshuffling (Metatibia)")

# 4. Protíbia na Floresta (~ Areia)
p4 <- ggplot(base_solo_for %>% drop_na(RS_Protibia, Areia), aes(x = Areia, y = RS_Protibia)) +
  geom_point(alpha = 0.7, size = 3, color = "#D55E00") +
  geom_smooth(method = "lm", color = "black", fill = "gray80") +
  theme_classic() +
  labs(title = "D) Forest Protibia (Soil Texture)",
       x = "Soil Sand Content (%)",
       y = "Functional Reshuffling (Protíbia)")

# Combinando os 4 gráficos em um único painel usando o pacote patchwork
painel_graficos <- (p1 | p2) / (p3 | p4)

# Exportando o painel em alta resolução (300 dpi) no padrão para artigos/teses
ggsave(here::here("Graficos_Modelos_Significativos.png"), plot = painel_graficos, 
       width = 11, height = 9, dpi = 600)

message("Gráficos exportados com sucesso! Verifique o arquivo 'Graficos_Modelos_Significativos.png'.")
