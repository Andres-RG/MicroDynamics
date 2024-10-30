# Determinación de indices de diversidad alfa raton 4 === --- ===
##---------------------------
# Se define la ruta de trabajo
setwd("MiceProject/")
# Librerias necesarias
library(vegan)
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/basal_period_mouse_4.RData")      # basal_subject4_aggregate
load(file = "03_out/data/fatdiet_period_mouse_4.RData")    # fatdiet_subject4_aggregate
load(file = "03_out/data/recovered1_period_mouse_4.RData") # recover1_subject4_aggregate
load(file = "03_out/data/vancomycin_period_mouse_4.RData") # vancomycin_subject4_aggregate
load(file = "03_out/data/recovered2_period_mouse_4.RData") # recover2_subject4_aggregate
load(file = "03_out/data/gentamicin_period_mouse_4.RData") # gentamicin_subject4_aggregate
load(file = "03_out/data/recovered3_period_mouse_4.RData") # recover3_subject4_aggregate
load(file = "03_out/data/mouse4.RData")                    # mouse4
##---------------------------
# Se guarda todas las observaciones del raton 2
mouse4 <- rbind(basal_subject4_aggregate,
                fatdiet_subject4_aggregate,
                recover1_subject4_aggregate,
                vancomycin_subject4_aggregate,
                recover2_subject4_aggregate,
                gentamicin_subject4_aggregate,
                recover3_subject4_aggregate)
# save(mouse4, file = "03_out/data/mouse4.RData")
##---------------------------
# Indices de diversidad
# simpson
simpson_mouse4 <- apply(mouse4, 1, function(x) diversity(x, index = "simpson"))
names(simpson_mouse4) <- c()
simpson_index_mouse4 <- data.frame(
  time = seq(1,length(simpson_mouse4),1),
  simpson = simpson_mouse4
)
rownames(simpson_index_mouse4) <- seq(1,length(simpson_mouse4),1)
simpson_index_mouse4
# save(simpson_index_mouse4, file = "03_out/data/index_diversity_simpson_mouse4.RData")
# berger-parker index
berger_parker_mouse4 <- apply(mouse4, 1,
                              function(x) max(x)/sum(x))
names(berger_parker_mouse4) <- c()
berger_parker_index_mouse4 <- data.frame(
  time = seq(1,length(berger_parker_mouse4),1),
  b_p = berger_parker_mouse4
)
rownames(berger_parker_index_mouse4) <- seq(1,length(berger_parker_mouse4),1)
berger_parker_index_mouse4
# save(berger_parker_index_mouse4, file = "03_out/data/index_diversity_berger_parker_mouse4.RData")
# shannon normalizado
shannon_mouse4 <- apply(mouse4,
                        1,
                        function(x) diversity(x, index = "shannon"))
S <- specnumber(mouse4)
shannon_normalized_mouse4 <- shannon_mouse4 / log(S)
names(shannon_normalized_mouse4) <- c()
shannon_normalized_index_mouse4 <- data.frame(
  time = seq(1,length(shannon_normalized_mouse4),1),
  shannon = shannon_normalized_mouse4
)
rownames(shannon_normalized_index_mouse4) <- seq(1,length(shannon_normalized_mouse4),1)
shannon_normalized_index_mouse4
# save(shannon_normalized_index_mouse4, file = "03_out/data/index_diversity_shannon_normalized_mouse4.RData")
# gini-simpson
gini_simpson_index_mouse4 <- data.frame(
  time = seq(1,length(simpson_mouse4),1),
  gini_simpson = 1 - simpson_index_mouse4[,2]
)
rownames(gini_simpson_index_mouse4) <- seq(1,length(simpson_mouse4),1)
gini_simpson_index_mouse4
# save(gini_simpson_index_mouse4, file = "03_out/data/index_diversity_gini_simpson_mouse4.RData")