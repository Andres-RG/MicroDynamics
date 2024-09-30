# Determinación de diversidad
library(vegan)
#
load(file = "03_out/data/basal_period_mouse_2.RData")      # basal_subject2_aggregate
load(file = "03_out/data/fatdiet_period_mouse_2.RData")    # fatdiet_subject2_aggregate
load(file = "03_out/data/recovered1_period_mouse_2.RData") # recover1_subject2_aggregate
load(file = "03_out/data/vancomycin_period_mouse_2.RData") # vancomycin_subject2_aggregate
load(file = "03_out/data/recovered2_period_mouse_2.RData") # recover2_subject2_aggregate
load(file = "03_out/data/gentamicin_period_mouse_2.RData") # gentamicin_subject2_aggregate
load(file = "03_out/data/recovered3_period_mouse_2.RData") # recover3_subject2_aggregate
mouse2 <- rbind(basal_subject2_aggregate,
                fatdiet_subject2_aggregate,
                recover1_subject2_aggregate,
                vancomycin_subject2_aggregate,
                recover2_subject2_aggregate,
                gentamicin_subject2_aggregate,
                recover3_subject2_aggregate)
# simpson
s_i <- apply(mouse2, 1, function(x) diversity(x, index = "simpson"))
names(s_i) <- c()
simpson_index <- data.frame(
  time = seq(1,length(s_i),1),
  simpson = s_i
)
rownames(simpson_index) <- seq(1,length(s_i),1)
simpson_index <- as.table(simpson_index)
simpson_index
# pielou
pielou_index <- readRDS("~/Documents/maestria/MicroDynamics/MiceProject/01_raw_data/mice_data.RDS")
pielou_index
# berger-parker index
b_p_index <- apply(mouse2, 1,
                             function(x) max(x)/sum(x))
names(b_p_index) <- c()
berger_parker_index <- data.frame(
  time = seq(1,length(b_p_index),1),
  b_p = b_p_index
)
rownames(berger_parker_index) <- seq(1,length(b_p_index),1)
berger_parker_index <- as.table(berger_parker_index)
berger_parker_index
