# ==============================================================================
# ESTRATÉGIA DE SIMPLIFICAÇÃO: MODELOS ADITIVOS (SEM INTERAÇÃO)
# ==============================================================================

# 1. Definir o Controle Local (PCA Ambiental Único)
controles_locais <- "+ PC1_Ambiente_Local"

lista_selecao_simples <- list()
lista_modelos_vencedores_simples <- list()

for(y in cols_metricas) {
  
  # Filtro de dados completos
  dados_modelagem <- matriz_analise %>% 
    filter(!is.na(!!sym(y)), !is.na(PC1_Ambiente_Local))
  
  modelos_candidatos <- list()
  
  # --- MODELO NULO ---
  modelos_candidatos[["Nulo"]] <- glm(as.formula(paste(y, "~ 1")), 
                                      data = dados_modelagem, family = gaussian)
  
  # --- LOOP DAS HIPÓTESES DE PAISAGEM ---
  for(nome_h in names(hipoteses_paisagem)) {
    var_x <- hipoteses_paisagem[[nome_h]]
    
    # VERSÃO A: Aditiva com Uso (Testa se a paisagem importa controlando pelo tipo de uso)
    # f_uso <- as.formula(paste(y, "~", var_x, "+ uso", controles_locais))
    
    # VERSÃO B: Paisagem Pura (A mais simples: Testa se a paisagem importa independente do uso)
    f_pura <- as.formula(paste(y, "~", var_x, controles_locais))
    
    try({
      # Vamos usar a f_pura primeiro, pois é a que mais economiza parâmetros (K baixo)
      mod <- glm(f_pura, data = dados_modelagem, family = gaussian, na.action = na.fail)
      modelos_candidatos[[nome_h]] <- mod
    }, silent = TRUE)
  }
  
  # Seleção de Modelos
  if(length(modelos_candidatos) > 0) {
    selecao <- model.sel(modelos_candidatos)
    df_sel <- as.data.frame(selecao)
    df_sel$Variavel_Resposta <- y
    df_sel$Hipoteses <- rownames(df_sel)
    lista_selecao_simples[[y]] <- df_sel
    
    melhor_nome <- rownames(df_sel)[1]
    lista_modelos_vencedores_simples[[y]] <- modelos_candidatos[[melhor_nome]]
  }
}

## ==============================================================================
# 2. EXPORTAR RESULTADOS SIMPLIFICADOS (CORRIGIDO)
# ==============================================================================

tabela_aic_simples <- bind_rows(lista_selecao_simples) %>%
  # Trocamos 'K' por 'df' e aproveitamos para renomear para 'K' se desejar
  rename(K = df) %>% 
  select(Variavel_Resposta, Hipoteses, K, AICc, delta, weight) %>%
  mutate(across(where(is.numeric), round, 3)) # Arredonda para facilitar a leitura

write.xlsx(tabela_aic_simples, "Resultado_AIC_Simplificado_Paisagem.xlsx")

# --- VERIFICAÇÃO RÁPIDA ---
# Olhe no console se algum modelo de paisagem superou o nulo (Delta > 2)
tabela_aic_simples %>% 
  filter(Variavel_Resposta == "SAD_Slope") %>% # Exemplo com uma métrica
  head()
