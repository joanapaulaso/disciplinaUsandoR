# CONFIGURANDO AMBIENTE DE TRABALHO
setwd("D:/aulasR_PPGBIO/projeto_final")
getwd()

# CARREGANDO OS PACOTES
library("reshape2")
library("trafo")
library("mdatools")
library("car")
library("vegan")
library("dplyr")
library("QuantPsyc")
library("ggplot2")
library("AMR")

# IMPORTANDO OS DADOS
all_data <- read.table("./r_data.txt", sep = "\t", header = TRUE)
all_data

# DEFININDO COLUNAS QUALITATIVAS E QUANTITATIVAS
desc_cols <- all_data[, 1:3]
rt_mz_cols <- all_data[, 2:3]
quant_cols <- all_data[, 4:26]

# TRANSPONDO A TABELA DE DADOS PARA EXECUTAR TESTES DE HOMOCEDASTICIDADE (Levene's) E NORMALIDADE DE DISTRIBUIÇÃO (Kolmogorov-Smirnov)
norm_data_melt <- melt(all_data, id.vars = c("ID", "Rt", "Mz"), measure.vars = c("QC1.NEG", "QC2.NEG", "QC3.NEG", "QC4.NEG", "QC5.NEG", "VA1EE1.NEG", "VA1EE2.NEG", "VA1EE3.NEG", "VA1EM1.NEG", "VA1EM2.NEG", "VA1EM3.NEG", "VB1EE1.NEG", "VB1EE2.NEG", "VB1EE3.NEG", "VB1EM1.NEG", "VB1EM2.NEG", "VB1EM3.NEG", "VC1EE1.NEG", "VC1EE2.NEG", "VC1EE3.NEG", "VC1EM1.NEG", "VC1EM2.NEG", "VC1EM3.NEG"))
leveneTest(norm_data_melt$value ~ norm_data_melt$variable, center = mean)
ks.test(norm_data_melt$value , y = "pnorm")

# EXTRAINDO AS COLUNAS DE CADA TIPO AMOSTRAL PARA NORMALIZAÇÃO
norm_data <- all_data[, 4:26]
qc <- all_data[, 4:8]
vaee <- all_data[, 9:11]
vaem <- all_data[, 12:14]
vbee <- all_data[, 15:17]
vbem <- all_data[, 18:20]
vcee <- all_data[, 21:23]
vcem <- all_data[, 24:26]

norm_data
norm_data_logss <- as.matrix(norm_data)

# TRANSFORMAÇÃO DE DADOS POR LOG10
vanilla_log <- log10(norm_data_logss)
plot(vanilla_log)
qqnorm(vanilla_log)
qqline(vanilla_log)
sum(is.infinite(vanilla_log))

# TRANSFORMAÇÃO DE DADOS POR GLOG
vanilla_lm <- lm(norm_data_melt$value ~ norm_data_melt$variable)
vanilla_glog <- glog(vanilla_lm)
glog_values <- vanilla_glog[["yt"]]
glog_values_df <- as.data.frame(glog_values)
norm_data_melt <- data.frame(c(norm_data_melt, glog_values_df))
norm_data_melt_glog <- norm_data_melt[, 1:4]
norm_data_melt_glog <- cbind(norm_data_melt_glog, norm_data_melt$glog_values)
colnames(norm_data_melt_glog)[5] <- "value"

# AUTOESCALONAMENTO DOS DADOS (NOVA TRANSPOSIÇÃO DA TABELA)
unmelted <- dcast(norm_data_melt_glog, ID ~ variable)
vanilla_autoscaled <- prep.autoscale(unmelted[, 2:ncol(unmelted)], center = TRUE, scale = TRUE)
vanilla_autoscaled <- cbind.data.frame(vanilla_autoscaled, desc_cols)
quant_vanilla_autoscaled <- na.omit(subset(vanilla_autoscaled, vanilla_autoscaled[, 1:23] >= 1))
quant_vanilla_autoscaled <- as.data.frame(quant_vanilla_autoscaled)
sum(is.na(quant_vanilla_autoscaled))

qc_norm <- quant_vanilla_autoscaled[, 1:5]
vaee_norm <- quant_vanilla_autoscaled[, 6:8]
vaem_norm <- quant_vanilla_autoscaled[, 9:11]
vbee_norm <- quant_vanilla_autoscaled[, 12:14]
vbem_norm <- quant_vanilla_autoscaled[, 15:17]
vcee_norm <- quant_vanilla_autoscaled[, 18:20]
vcem_norm <- quant_vanilla_autoscaled[, 21:23]

# CÁLCULO DA MÉDIA DAS REPLICATAS NORMALIZADAS
quant_vanilla_autoscaled$qc <- rowMeans(qc_norm)
quant_vanilla_autoscaled$vaee <- rowMeans(vaee_norm)
quant_vanilla_autoscaled$vaem <- rowMeans(vaem_norm)
quant_vanilla_autoscaled$vbee <- rowMeans(vbee_norm)
quant_vanilla_autoscaled$vbem <- rowMeans(vbem_norm)
quant_vanilla_autoscaled$vcee <- rowMeans(vcee_norm)
quant_vanilla_autoscaled$vcem <- rowMeans(vcem_norm)

# CÁLCULO DO DESVIO PADRÃO DAS REPLICARAS NORMALIZADAS
qc_sd <- apply(qc_norm, 1, sd)
vaee_sd <- apply(vaee_norm, 1, sd)
vaem_sd <- apply(vaem_norm, 1, sd)
vbee_sd <- apply(vbee_norm, 1, sd)
vbem_sd <- apply(vbem_norm, 1, sd)
vcee_sd <- apply(vcee_norm, 1, sd)
vcem_sd <- apply(vcem_norm, 1, sd)

quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, qc_sd = qc_sd)
quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, vaee_sd = vaee_sd)
quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, vaem_sd = vaem_sd)
quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, vbee_sd = vbee_sd)
quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, vbem_sd = vbem_sd)
quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, vcee_sd = vcee_sd)
quant_vanilla_autoscaled <- cbind(quant_vanilla_autoscaled, vcem_sd = vcem_sd)

# CÁLCULO DO COEFICIENTE DE VARIAÇÃO DAS REPLICATAS NORMALIZADAS
quant_vanilla_autoscaled$qc_cv <- quant_vanilla_autoscaled$qc_sd / quant_vanilla_autoscaled$qc
quant_vanilla_autoscaled$vaee_cv <- quant_vanilla_autoscaled$vaee_sd / quant_vanilla_autoscaled$vaee
quant_vanilla_autoscaled$vaem_cv <- quant_vanilla_autoscaled$vaem_sd / quant_vanilla_autoscaled$vaem
quant_vanilla_autoscaled$vbee_cv <- quant_vanilla_autoscaled$vbee_sd / quant_vanilla_autoscaled$vbee
quant_vanilla_autoscaled$vbem_cv <- quant_vanilla_autoscaled$vbem_sd / quant_vanilla_autoscaled$vbem
quant_vanilla_autoscaled$vcee_cv <- quant_vanilla_autoscaled$vcee_sd / quant_vanilla_autoscaled$vcee
quant_vanilla_autoscaled$vcem_cv <- quant_vanilla_autoscaled$vcem_sd / quant_vanilla_autoscaled$vcem

names(quant_vanilla_autoscaled)

# ANÁLISE DE COMPONENTES PRINCIPAIS
vanilla_pca <- prcomp(quant_vanilla_autoscaled[, 27:33])
biplot(vanilla_pca)
ggplot_pca(vanilla_pca, base_textsize = 15, arrows_colour = "#ad0252", arrows_textsize = 5)
bstick(vanilla_pca)
eigenvals(vanilla_pca)
summary(vanilla_pca)
loadings(vanilla_pca)
scores(vanilla_pca)

# EXTRAÇÃO DOS DADOS PARA EXECUÇÃO DA RDA
va <- grepl("va", colnames(quant_vanilla_autoscaled))
va <- quant_vanilla_autoscaled[va]
va$avg <- rowMeans(va)
va <- va$avg

vb <- grepl("vb", colnames(quant_vanilla_autoscaled))
vb <- quant_vanilla_autoscaled[vb]
vb$avg <- rowMeans(vb)
vb <- vb$avg

vc <- grepl("vc", colnames(quant_vanilla_autoscaled))
vc <- quant_vanilla_autoscaled[vc]
vc$avg <- rowMeans(vc)
vc <- vc$avg

ee <- grepl("ee", colnames(quant_vanilla_autoscaled))
ee <- quant_vanilla_autoscaled[ee]
ee$avg <- rowMeans(ee)
ee <- ee$avg

em <- grepl("em", colnames(quant_vanilla_autoscaled))
em <- quant_vanilla_autoscaled[em]
em$avg <- rowMeans(em)
em <- em$avg

v_avg_samples <- grepl("v", colnames(quant_vanilla_autoscaled))
v_avg_samples <- quant_vanilla_autoscaled[v_avg_samples]
v_avg_samples <- v_avg_samples[, 1:6]

rda_data <- cbind.data.frame(va, vb, vc, ee, em)
rda_data_species <- cbind.data.frame(va, vb, vc)
rda_data_treatments <- cbind.data.frame(ee, em)

# ANÁLISE DE REDUNDÂNCIA (RDA)
v_rda_t <- rda(v_avg_samples ~ ., data = rda_data_treatments)
ordiplot(v_rda_t, scaling = 2, main = "RDA - effect of extraction on samples")
summary(v_rda_t)

v_rda_s <- rda(v_avg_samples ~ ., data = rda_data_species)
ordiplot(v_rda_s, scaling = 2, main = "RDA - effect of species on samples")
summary(v_rda_s)

v_rda_a <- rda(v_avg_samples ~ ., data = rda_data)
ordiplot(v_rda_a, scaling = 2, main = "RDA - effect all groups on samples")
summary(v_rda_a)

# SELEÇÃO DE VARIÁVEIS SIGNIFICATIVAS DA RDA
rda_sel_t <- ordiR2step(rda(v_avg_samples ~ 1, data = rda_data_treatments), 
                      scope = formula(v_rda_t), 
                      R2scope = TRUE, 
                      steps = 1000,
                      trace = TRUE) 

rda_sel_t$call
selected_t <- rda(formula = v_avg_samples ~ ee + em, data = rda_data_treatments)
RsquareAdj(selected_t)
anova.cca(selected_t, step = 1000, by = "term")
ordiplot(selected_t, scaling = 2, type = "text")

rda_sel_s <- ordiR2step(rda(v_avg_samples ~ 1, data = rda_data_species),
                      scope = formula(v_rda_s),
                      R2scope = TRUE,
                      steps = 1000,
                      trace = TRUE)

rda_sel_s$call
selected_s <- rda(formula = v_avg_samples ~ va + vc + vb, data = rda_data_species)
RsquareAdj(selected_s)
anova.cca(selected_s, step = 1000, by = "term")
ordiplot(selected_s, scaling = 2, type = "text")

rda_sel_a <- ordiR2step(rda(v_avg_samples ~ 1, data = rda_data),
                      scope = formula(v_rda_a),
                      R2scope = TRUE,
                      steps = 1000,
                      trace = TRUE)

rda_sel_a$call
selected_a <- rda(formula = v_avg_samples ~ va + vc + vb + ee, data = rda_data)
RsquareAdj(selected_a)
anova.cca(selected_a, step = 1000, by = "term")
ordiplot(selected_a, scaling = 2, type = "text")

# TRANSPOSIÇÃO DOS DADOS PARA OS TESTES DE NORMALIDADE DE DISTRIUIÇÃO E HOMOCEDASTICIDADE
v_autoscaled_melt <- melt(vanilla_autoscaled, id.vars = c(24, 25, 26), measure.vars = c("QC1.NEG", "QC2.NEG", "QC3.NEG", "QC4.NEG", "QC5.NEG", "VA1EE1.NEG", "VA1EE2.NEG", "VA1EE3.NEG", "VA1EM1.NEG", "VA1EM2.NEG", "VA1EM3.NEG", "VB1EE1.NEG", "VB1EE2.NEG", "VB1EE3.NEG", "VB1EM1.NEG", "VB1EM2.NEG", "VB1EM3.NEG", "VC1EE1.NEG", "VC1EE2.NEG", "VC1EE3.NEG", "VC1EM1.NEG", "VC1EM2.NEG", "VC1EM3.NEG"))

norm_ss <- subset(v_autoscaled_melt, v_autoscaled_melt$value >= -1)
norm_ss$std <- as.numeric(decostand(norm_ss$value, "max"))
norm_ss_noqc <- subset(norm_ss, grepl("V", norm_ss$variable, fixed = TRUE))

sp_groups <- ifelse(grepl("VA", norm_ss_noqc$variable, fixed = TRUE), "VA", 
                    ifelse(grepl("VB", norm_ss_noqc$variable, fixed = TRUE), "VB", "VC"))
sp_groups <- as.data.frame(sp_groups)

treat_groups <- ifelse(grepl("EE", norm_ss_noqc$variable, fixed = TRUE), "EE", 
                    ifelse(grepl("EM", norm_ss_noqc$variable, fixed = TRUE), "EM", NA))
treat_groups <- as.data.frame(treat_groups)

norm_ss_noqc_g <- na.omit(cbind.data.frame(norm_ss_noqc, sp_groups, treat_groups))

hist(v_autoscaled_melt$value)
hist(norm_ss_noqc_g$std)
qqnorm(v_autoscaled_melt$value)
qqline(v_autoscaled_melt$value)
qqnorm(norm_ss_noqc_g$std)
qqline(norm_ss_noqc_g$std)
leveneTest(norm_ss_noqc_g$std ~ norm_ss_noqc_g$variable)
ks.test(norm_ss_noqc_g$std , y = "pnorm")
plot(norm_ss_noqc_g$value ~ norm_ss_noqc_g$variable)

param_aov_sp <- aov(norm_ss_noqc_g$value ~ norm_ss_noqc_g$sp_groups)
summary(param_aov_sp)
param_tukey_sp <- TukeyHSD(param_aov_sp)
param_tukey_sp
plot(param_tukey_sp)

param_aov_treat <- aov(norm_ss_noqc_g$value ~ norm_ss_noqc_g$treat_groups)
summary(param_aov_treat)
param_tukey_treat <- TukeyHSD(param_aov_treat)
param_tukey_treat
plot(param_tukey_treat)

kruskal.test(norm_ss_noqc_g$value ~ norm_ss_noqc_g$sp_groups)
pairwise.wilcox.test(norm_ss_noqc_g$value, norm_ss_noqc_g$sp_groups, p.adjust.method = "BH")

kruskal.test(norm_ss_noqc_g$value ~ norm_ss_noqc_g$treat_groups)
pairwise.wilcox.test(norm_ss_noqc_g$value, norm_ss_noqc_g$treat_groups, p.adjust.method = "BH")
t.test(norm_ss_noqc_g$value ~ norm_ss_noqc_g$treat_groups)

# TESTANDO A DISTRIBUIÇÃO, HOMOCEDASTICIDADE DOS SINAIS NOS GRUPOS DE ESPÉCIES
id_sp <- subset(norm_ss_noqc_g, norm_ss_noqc_g$ID == 11624)
leveneTest(id_sp$value ~ id_sp$sp_groups)
shapiro.test(id_sp$value)

# TESTANDO A DISTRIBUIÇÃO, HOMOCEDASTICIDADE DOS SINAIS NOS GRUPOS DE EXTRAÇÕES
id_treat <- subset(norm_ss_noqc_g, norm_ss_noqc_g$ID == 1790)
leveneTest(id_treat$value ~ id_treat$sp_groups)
shapiro.test(id_treat$value)

# REALIZANDO ANOVA E TUKEY NO CASO DE APROVADOS
summary(aov(id_sp$value ~ id_sp$sp_groups))
idTukey_sp <- TukeyHSD(aov(id_sp$value ~ id_sp$sp_groups))
idTukey_sp
plot(idTukey_sp)

summary(aov(id_treat$value ~ id_treat$sp_groups))
idTukey_treat <- TukeyHSD(aov(id_treat$value ~ id_treat$sp_groups))
idTukey_treat
plot(idTukey_treat)

# REALIZANDO KRUSKAL-WALLIS E WILLCOX NO CASO DE REJEITADOS
kruskal.test(id_sp$value ~ id_sp$sp_groups)
pairwise.wilcox.test(id_sp$value ~ id_sp$sp_groups, p.adjust.method = "BH")
boxplot(id_sp$value ~ id_sp$sp_groups)

kruskal.test(id_treat$value ~ id_treat$sp_groups)
pairwise.wilcox.test(id_sp$value ~ id_treat$sp_groups, p.adjust.method = "BH")

# SALVANDO OS RESULTADOS
write.csv(norm_ss_noqc_g, "var.csv", row.names = TRUE)
