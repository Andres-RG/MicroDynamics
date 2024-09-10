# RATON 2 === --- === subject2 (en los datos)
##---------------------------
# Librerias necesarias
library(readr)
library(mlBioNets)
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/asv_table.RData") # asv_table
load(file = "03_out/data/asv_table_aggregate.RData") # asv_table_cllps
load(file = "03_out/data/basal_period_mouse_2.RData")      # basal_subject2_aggregate
load(file = "03_out/data/fatdiet_period_mouse_2.RData")    # fatdiet_subject2_aggregate
load(file = "03_out/data/recovered1_period_mouse_2.RData") # recover1_subject2_aggregate
load(file = "03_out/data/vancomycin_period_mouse_2.RData") # vancomycin_subject2_aggregate
load(file = "03_out/data/recovered2_period_mouse_2.RData") # recover2_subject2_aggregate
load(file = "03_out/data/gentamicin_period_mouse_2.RData") # gentamicin_subject2_aggregate
load(file = "03_out/data/recovered3_period_mouse_2.RData") # recover3_subject2_aggregate
##---------------------------
# Se cargan los datos del script
source("02_scripts/01_procesamiento_y_filtrado_de_datos.R")
##---------------------------
# De cada raton, se separan las abundancias por tratamiento
## RATON 2 ##----------------
## Basal - --- --- --- --- --
basal_start_t            <- subject2$time[1]
fatdiet_start_t          <- asv_pert[which(asv_pert$subject == 2),]$start[1]
basal_end_t              <- subject2$time[which(subject2$time == fatdiet_start_t) -1 ]
basal_subject2_IDs       <- subject2[which(subject2$time == basal_start_t) : which(subject2$time == basal_end_t),]$sampleID
basal_subject2_aggregate <- asv_table_aggregate[basal_subject2_IDs,]
#   Se obtienen los datos agregados del periodo basal del raton 2
#   Se guarda el objeto procesado como .RData
# save(basal_subject2_aggregate, file = "03_out/data/basal_period_mouse_2.RData")
## Fat diet --- --- --- --- -
fatdiet_start_t
fatdiet_end_t              <- asv_pert[which(asv_pert$subject == 2),]$end[1]
fatdiet_subject2_IDs       <- subject2[which(subject2$time == fatdiet_start_t) : which(subject2$time == fatdiet_end_t),]$sampleID
fatdiet_subject2_aggregate <- asv_table_aggregate[fatdiet_subject2_IDs,]
#   Se obtienen los datos agregados del periodo fatdiet del raton 2
#   Se guarda el objeto procesado como .RData
# save(fatdiet_subject2_aggregate, file = "03_out/data/fatdiet_period_mouse_2.RData")
## Recovered 1 -- --- --- ---
recover1_start_t            <- subject2$time[which(subject2$time == fatdiet_end_t)+1]
vancomycin_start_t          <- asv_pert[which(asv_pert$subject == 2),]$start[2]
recover1_end_t              <- subject2$time[which(subject2$time == vancomycin_start_t)-1]
recover1_subject2_IDs       <- subject2[which(subject2$time == recover1_start_t) : which(subject2$time == recover1_end_t),]$sampleID
recover1_subject2_aggregate <- asv_table_aggregate[recover1_subject2_IDs,]
#   Se obtienen los datos agregados del periodo de recuperacion 1 del raton 2
#   Se guarda el objeto procesado como .RData
# save(recover1_subject2_aggregate, file = "03_out/data/recovered1_period_mouse_2.RData")
## Vancomycin --- --- --- ---
vancomycin_start_t
vancomycin_end_t              <- asv_pert[which(asv_pert$subject == 2),]$end[2]
vancomycin_subject2_IDs       <- subject2[which(subject2$time == vancomycin_start_t) : which(subject2$time == vancomycin_end_t),]$sampleID
vancomycin_subject2_aggregate <- asv_table_aggregate[vancomycin_subject2_IDs,]
#   Se obtienen los datos agregados del periodo con vancomicina del raton 2
#   Se guarda el objeto procesado como .RData
# save(vancomycin_subject2_aggregate, file = "03_out/data/vancomycin_period_mouse_2.RData")
## Recovered 2 -- --- --- ---
recover2_start_t            <- subject2$time[which(subject2$time == vancomycin_end_t)+1]
gentamicin_start_t          <- asv_pert[which(asv_pert$subject == 2),]$start[3]
recover2_end_t              <- subject2$time[which(subject2$time == gentamicin_start_t)-1]
recover2_subject2_IDs       <- subject2[which(subject2$time == recover2_start_t) : which(subject2$time == recover2_end_t),]$sampleID
recover2_subject2_aggregate <- asv_table_aggregate[recover2_subject2_IDs,]
#   Se obtienen los datos agregados del periodo de recuperacion 2 del raton 2
#   Se guarda el objeto procesado como .RData
# save(recover2_subject2_aggregate, file = "03_out/data/recovered2_period_mouse_2.RData")
gentamicin_start_t
gentamicin_end_t              <- asv_pert[which(asv_pert$subject == 2),]$end[3]
gentamicin_subject2_IDs       <- subject2[which(subject2$time == gentamicin_start_t) : which(subject2$time == gentamicin_end_t),]$sampleID
gentamicin_subject2_aggregate <- asv_table_aggregate[gentamicin_subject2_IDs,]
#   Se obtienen los datos agregados del periodo con gentamicina del raton 2
#   Se guarda el objeto procesado como .RData
# save(gentamicin_subject2_aggregate, file = "03_out/data/gentamicin_period_mouse_2.RData")
## Recovered 3 -- --- --- ---
recover3_start_t            <- subject2$time[which(subject2$time == gentamicin_end_t)+1]
recover3_end_t              <- subject2$time[which.max(subject2$time)]
recover3_subject2_IDs       <- subject2[which(subject2$time == recover3_start_t) : which(subject2$time == recover3_end_t),]$sampleID
recover3_subject2_aggregate <- asv_table_aggregate[recover3_subject2_IDs,]
#   Se obtienen los datos agregados del periodo de recuperacion 3 del raton 2
#   Se guarda el objeto procesado como .RData
# save(recover3_subject2_aggregate, file = "03_out/data/recovered3_period_mouse_2.RData")
