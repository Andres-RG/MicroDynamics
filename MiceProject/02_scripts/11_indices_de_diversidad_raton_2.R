# Determinación de indices de diversidad alfa raton 2 === --- ===
##---------------------------
# Se define la ruta de trabajo
setwd("MiceProject/")
# Librerias necesarias
library(vegan)
  ##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/basal_period_mouse_2.RData")      # basal_subject2_aggregate
load(file = "03_out/data/fatdiet_period_mouse_2.RData")    # fatdiet_subject2_aggregate
load(file = "03_out/data/recovered1_period_mouse_2.RData") # recover1_subject2_aggregate
load(file = "03_out/data/vancomycin_period_mouse_2.RData") # vancomycin_subject2_aggregate
load(file = "03_out/data/recovered2_period_mouse_2.RData") # recover2_subject2_aggregate
load(file = "03_out/data/gentamicin_period_mouse_2.RData") # gentamicin_subject2_aggregate
load(file = "03_out/data/recovered3_period_mouse_2.RData") # recover3_subject2_aggregate
load(file = "03_out/data/mouse2.RData")                    # mouse2
##---------------------------
# Se guarda todas las observaciones del raton 2
mouse2 <- rbind(basal_subject2_aggregate,
                fatdiet_subject2_aggregate,
                recover1_subject2_aggregate,
                vancomycin_subject2_aggregate,
                recover2_subject2_aggregate,
                gentamicin_subject2_aggregate,
                recover3_subject2_aggregate)
# save(mouse2, file = "03_out/data/mouse2.RData")
##---------------------------
# Indices de diversidad
# simpson
# probabilidad de que dos individuos de una comunidad seleccionados al azar,
# pertenezcan a la misma especie. Valores altos = menor diversidad
simpson_mouse2 <- apply(mouse2, 1, function(x) diversity(x, index = "simpson"))
names(simpson_mouse2) <- c()
simpson_index_mouse2 <- data.frame(
  time = seq(1,length(simpson_mouse2),1),
  simpson = simpson_mouse2
)
rownames(simpson_index_mouse2) <- seq(1,length(simpson_mouse2),1)
simpson_index_mouse2
# save(simpson_index_mouse2, file = "03_out/data/index_diversity_simpson_mouse2.RData")
# pielou
pielou_index_mouse2 <- readRDS("~/Documents/maestria/MicroDynamics/MiceProject/01_raw_data/mice_data.RDS")
pielou_index_mouse2
# save(pielou_index_mouse2, file = "03_out/data/index_diversity_pielou_mouse2.RData")
# berger-parker index
berger_parker_mouse2 <- apply(mouse2, 1,
                             function(x) max(x)/sum(x))
names(berger_parker_mouse2) <- c()
berger_parker_index_mouse2 <- data.frame(
  time = seq(1,length(berger_parker_mouse2),1),
  b_p = berger_parker_mouse2
)
rownames(berger_parker_index_mouse2) <- seq(1,length(berger_parker_mouse2),1)
berger_parker_index_mouse2
# save(berger_parker_index_mouse2, file = "03_out/data/index_diversity_berger_parker_mouse2.RData")
# gini-simpson
# qué tan uniformemente se distribuyen los individuos en las diferentes especies
# de una comunidad. Valor alto = están igualmente representadas. Valor bajo = algunas especies dominan
gini_simpson_index_mouse2 <- data.frame(
  time = seq(1,length(simpson_mouse2),1),
  gini_simpson = 1 - simpson_index_mouse2[,2]
)
rownames(gini_simpson_index_mouse2) <- seq(1,length(simpson_mouse2),1)
gini_simpson_index_mouse2
# save(gini_simpson_index_mouse2, file = "03_out/data/index_diversity_gini_simpson_mouse2.RData")
