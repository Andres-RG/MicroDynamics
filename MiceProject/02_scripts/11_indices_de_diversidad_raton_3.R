# Determinación de indices de diversidad alfa raton 3 === --- ===
##---------------------------
# Se define la ruta de trabajo
setwd("MiceProject/")
# Librerias necesarias
library(vegan)
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/basal_period_mouse_3.RData")      # basal_subject3_aggregate
load(file = "03_out/data/fatdiet_period_mouse_3.RData")    # fatdiet_subject3_aggregate
load(file = "03_out/data/recovered1_period_mouse_3.RData") # recover1_subject3_aggregate
load(file = "03_out/data/vancomycin_period_mouse_3.RData") # vancomycin_subject3_aggregate
load(file = "03_out/data/recovered2_period_mouse_3.RData") # recover2_subject3_aggregate
load(file = "03_out/data/gentamicin_period_mouse_3.RData") # gentamicin_subject3_aggregate
load(file = "03_out/data/recovered3_period_mouse_3.RData") # recover3_subject3_aggregate
load(file = "03_out/data/mouse3.RData")                    # mouse3
##---------------------------
# Se guarda todas las observaciones del raton 2
mouse3 <- rbind(basal_subject3_aggregate,
                fatdiet_subject3_aggregate,
                recover1_subject3_aggregate,
                vancomycin_subject3_aggregate,
                recover2_subject3_aggregate,
                gentamicin_subject3_aggregate,
                recover3_subject3_aggregate)
# save(mouse3, file = "03_out/data/mouse3.RData")
##---------------------------
# Indices de diversidad
# simpson
simpson_mouse3 <- apply(mouse3, 1, function(x) diversity(x, index = "simpson"))
names(simpson_mouse3) <- c()
simpson_index_mouse3 <- data.frame(
  time = seq(1,length(simpson_mouse3),1),
  simpson = simpson_mouse3
)
rownames(simpson_index_mouse3) <- seq(1,length(simpson_mouse3),1)
simpson_index_mouse3
# save(simpson_index_mouse3, file = "03_out/data/index_diversity_simpson_mouse3.RData")
# berger-parker index
berger_parker_mouse3 <- apply(mouse3, 1,
                              function(x) max(x)/sum(x))
names(berger_parker_mouse3) <- c()
berger_parker_index_mouse3 <- data.frame(
  time = seq(1,length(berger_parker_mouse3),1),
  b_p = berger_parker_mouse3
)
rownames(berger_parker_index_mouse3) <- seq(1,length(berger_parker_mouse3),1)
berger_parker_index_mouse3
# save(berger_parker_index_mouse3, file = "03_out/data/index_diversity_berger_parker_mouse3.RData")
# shannon normalizado
shannon_mouse3 <- apply(mouse3,
                        1,
                        function(x) diversity(x, index = "shannon"))
S <- apply(mouse3, 1, function(x) sum(x > 0))
shannon_normalized_mouse3 <- shannon_mouse3 / log(S)
names(shannon_normalized_mouse3) <- c()
shannon_normalized_index_mouse3 <- data.frame(
  time = seq(1,length(shannon_normalized_mouse3),1),
  shannon = shannon_normalized_mouse3
)
rownames(shannon_normalized_index_mouse3) <- seq(1,length(shannon_normalized_mouse3),1)
shannon_normalized_index_mouse3
# save(shannon_normalized_index_mouse3, file = "03_out/data/index_diversity_shannon_normalized_mouse3.RData")
# pielou
shannon_mouse3
S
pielou_mouse3 <- shannon_mouse3 / S
names(pielou_mouse3) <-  c()
pielou_index_mouse3 <-  data.frame(
  time = seq(1,length(pielou_mouse3), 1),
  pielou = pielou_mouse3
)
pielou_index_mouse3
# save(pielou_index_mouse3, file = "03_out/data/index_diversity_pielou_mouse3.RData")
# gini-simpson
gini_simpson_index_mouse3 <- data.frame(
  time = seq(1,length(simpson_mouse3),1),
  gini_simpson = 1 - simpson_index_mouse3[,2]
)
rownames(gini_simpson_index_mouse3) <- seq(1,length(simpson_mouse3),1)
gini_simpson_index_mouse3
# save(gini_simpson_index_mouse3, file = "03_out/data/index_diversity_gini_simpson_mouse3.RData")
