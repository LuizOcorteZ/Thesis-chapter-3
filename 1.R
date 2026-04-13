# --- SUPLEMENTAR 1: CARACTERIZAÇÃO DO CONTROLE LOCAL ---
library(factoextra)
library(openxlsx)

# 1. Tabela de Loadings Detalhada
# Extrai os pesos das variáveis para os 3 primeiros eixos
loadings_detalhado <- as.data.frame(pca_amb$rotation[, 1:3])
write.xlsx(loadings_detalhado, "Tabela_S1_Loadings_Ambientais.xlsx", rowNames = TRUE)

# 2. Gráfico Biplot (O visual principal do suplementar)
# Mostra como temperatura e vegetação se opõem no espaço multivariado
p_biplot <- fviz_pca_var(pca_amb,
                         col.var = "contrib", 
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE,
                         title = "Figura S1: Gradiente de Estresse Térmico e Estrutura Local") +
  theme_minimal()

# 3. Gráfico de Variância (Scree Plot)
# Demonstra que o PC1 é o eixo dominante (provavelmente > 50% da variância)
p_scree <- fviz_eig(pca_amb, addlabels = TRUE, ylim = c(0, 80),
                    title = "Figura S2: Variância Explicada pelo PCA Ambiental")

# Salvar arquivos
ggsave("Figura_S1_Biplot_Ambiental.png", p_biplot, width = 8, height = 7)
ggsave("Figura_S2_ScreePlot.png", p_scree, width = 6, height = 5)

# --- SUPLEMENTAR 2: SELEÇÃO DE MODELOS DETALHADA ---
library(MuMIn)
library(openxlsx)

# 1. Função para consolidar as tabelas de seleção do loop anterior
# Assume que 'lista_selecao' contém os dataframes gerados pelo model.sel()
tabela_suplementar_modelos <- bind_rows(lista_selecao) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  select(Variavel_Resposta, Hipoteses, df, logLik, AICc, delta, weight) %>%
  rename(Delta_AICc = delta, Peso_AICc = weight)

# 2. Salvar em Excel
write.xlsx(tabela_suplementar_modelos, "Tabela_S2_Ranking_AIC_Completo.xlsx")

# 3. Gráfico de Resíduos dos Modelos Vencedores
# Importante para provar que as premissas estatísticas foram atendidas
pdf("Figura_S3_Diagnostico_Residuos.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
for(y in names(lista_modelos_vencedores)) {
  mod <- lista_modelos_vencedores[[y]]
  plot(mod, which = 1, main = paste("Resíduos vs Ajustados:", y))
  plot(mod, which = 2, main = paste("Normal Q-Q:", y))
}
dev.off()