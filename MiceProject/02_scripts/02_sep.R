# Sujeto 2
##-----------------
# Librerias necesarias
library(readr)
library(mlBioNets)
##-----------------
# Se cargan los datos procesados
load(file = "03_out/data/asv_table.RData") # asv_table
load(file = "03_out/data/asv_tabla_collapsed.RData") # asv_table_cllps
##-----------------
# se cargan los datos anteriores
source("02_scripts/01_procesamiento_y_filtrado_de_datos.R")
##-----------------
# Se seleccionan los datos únicamente del ratón 2. Es al que se le aplicó un
# tratamiento completo
subject2 <- asv_meta[which(asv_meta$subject == 2),]
S2<-S2[order(S2$time), ]