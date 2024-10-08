# Determinación de indices de diversidad alfa raton 5 === --- ===
##---------------------------
# Se define la ruta de trabajo
setwd("MiceProject/")
# Librerias necesarias
library(vegan)
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/basal_period_mouse_5.RData")      # basal_subject5_aggregate
load(file = "03_out/data/fatdiet_period_mouse_5.RData")    # fatdiet_subject5_aggregate
load(file = "03_out/data/recovered1_period_mouse_5.RData") # recover1_subject5_aggregate
load(file = "03_out/data/vancomycin_period_mouse_5.RData") # vancomycin_subject5_aggregate
load(file = "03_out/data/recovered2_period_mouse_5.RData") # recover2_subject5_aggregate
load(file = "03_out/data/gentamicin_period_mouse_5.RData") # gentamicin_subject5_aggregate
load(file = "03_out/data/recovered3_period_mouse_5.RData") # recover3_subject5_aggregate
load(file = "03_out/data/mouse5.RData")                    # mouse5
##---------------------------
# Se guarda todas las observaciones del raton 2
mouse5 <- rbind(basal_subject5_aggregate,
                fatdiet_subject5_aggregate,
                recover1_subject5_aggregate,
                vancomycin_subject5_aggregate,
                recover2_subject5_aggregate,
                gentamicin_subject5_aggregate,
                recover3_subject5_aggregate)
# save(mouse5, file = "03_out/data/mouse5.RData")
##---------------------------
# Indices de diversidad
# simpson
simpson_mouse5 <- apply(mouse5, 1, function(x) diversity(x, index = "simpson"))
names(simpson_mouse5) <- c()
simpson_index_mouse5 <- data.frame(
  time = seq(1,length(simpson_mouse5),1),
  simpson = simpson_mouse5
)
rownames(simpson_index_mouse5) <- seq(1,length(simpson_mouse5),1)
simpson_index_mouse5
# save(simpson_index_mouse5, file = "03_out/data/index_diversity_simpson_mouse5.RData")
# berger-parker index
berger_parker_mouse5 <- apply(mouse5, 1,
                              function(x) max(x)/sum(x))
names(berger_parker_mouse5) <- c()
berger_parker_index_mouse5 <- data.frame(
  time = seq(1,length(berger_parker_mouse5),1),
  b_p = berger_parker_mouse5
)
rownames(berger_parker_index_mouse5) <- seq(1,length(berger_parker_mouse5),1)
berger_parker_index_mouse5
# save(berger_parker_index_mouse5, file = "03_out/data/index_diversity_berger_parker_mouse5.RData")
# shannon normalizado
shannon_mouse5 <- apply(mouse5,
                        1,
                        function(x) diversity(x, index = "shannon"))
S <- apply(mouse5, 1, function(x) sum(x > 0))
shannon_normalized_mouse5 <- shannon_mouse5 / log(S)
names(shannon_normalized_mouse5) <- c()
shannon_normalized_index_mouse5 <- data.frame(
  time = seq(1,length(shannon_normalized_mouse5),1),
  shannon = shannon_normalized_mouse5
)
rownames(shannon_normalized_index_mouse5) <- seq(1,length(shannon_normalized_mouse5),1)
shannon_normalized_index_mouse5
# save(shannon_normalized_index_mouse5, file = "03_out/data/index_diversity_shannon_normalized_mouse5.RData")
# pielou
shannon_mouse5
S
pielou_mouse5 <- shannon_mouse5 / S
names(pielou_mouse5) <-  c()
pielou_index_mouse5 <-  data.frame(
  time = seq(1,length(pielou_mouse5), 1),
  pielou = pielou_mouse5
)
pielou_index_mouse5
# save(pielou_index_mouse5, file = "03_out/data/index_diversity_pielou_mouse5.RData")
# gini-simpson
gini_simpson_index_mouse5 <- data.frame(
  time = seq(1,length(simpson_mouse5),1),
  gini_simpson = 1 - simpson_index_mouse5[,2]
)
rownames(gini_simpson_index_mouse5) <- seq(1,length(simpson_mouse5),1)
gini_simpson_index_mouse5
# save(gini_simpson_index_mouse5, file = "03_out/data/index_diversity_gini_simpson_mouse5.RData")