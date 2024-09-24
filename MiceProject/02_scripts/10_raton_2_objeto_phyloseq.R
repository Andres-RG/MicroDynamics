# Phyloseq object mouse 2 === --- ===
##---------------------------
# Se define la ruta
setwd("MiceProject/")
##---------------------------
# Librerias necesarias
library(devtools)
library(readr)
library(mlBioNets)
library(igraph)
library(minet)
library(ggplot2)
library(gridExtra)
library(phyloseq)
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/asv_table.RData") # asv_table
load(file = "03_out/data/asv_table_aggregate.RData") # asv_table_cllps
##---------------------------
# Se cargan los datos procesados
source(file = "02_scripts/01_procesamiento_y_filtrado_de_datos.R")
##---------------------------
head(asv_taxa) # creo que solo es del raton 2
asv_taxa_simplificada <- asv_taxa[,-c(1,2)] # del asv_taxa, se eliminan las columnas 1 y 2, que tienen
                                            # el asv_number y la secuencia
rownames(asv_taxa_simplificada) <- asv_ids # se colocan los asv_ids como nombre en los renglones
head(asv_taxa_simplificada)
str(asv_taxa_simplificada)
otu_table_complete_s2 <- otu_table(asv_table, taxa_are_rows = T) # se genera la otu_table
# contiene los conteos de cada asv en cada tiempo
##---------------------------
sample_IDs_s2 <- c(rownames(basal_subject2_aggregate),
                rownames(fatdiet_subject2_aggregate),
                rownames(recover1_subject2_aggregate),
                rownames(vancomycin_subject2_aggregate),
                rownames(recover2_subject2_aggregate),
                rownames(gentamicin_subject2_aggregate),
                rownames(recover3_subject2_aggregate))
treatments_s2 <- c(rep("Basal", length(rownames(basal_subject2_aggregate))),
                rep("Fat diet", length(rownames(fatdiet_subject2_aggregate))),
                rep("Recovered 1", length(rownames(recover1_subject2_aggregate))),
                rep("Vancomycin", length(rownames(vancomycin_subject2_aggregate))),
                rep("Recovered 2", length(rownames(recover2_subject2_aggregate))),
                rep("Gentamycin", length(rownames(gentamicin_subject2_aggregate))),
                rep("Recovered 3", length(rownames(recover3_subject2_aggregate))))
##---------------------------
meta_samples <-  data.frame(
  "samples" = sample_IDs_s2,
  "treatments" = treatments_s2
)
row.names(meta_samples) <- sample_IDs_s2
# se genera una base de datos de los datos de las muestras, que tienen los id's de las muestras
# y el periodo de experimento al que pertenecen
# creo que solo es del raton 2
##---------------------------
otu_table_complete_s2 <- otu_table_complete_s2[,which(colnames(otu_table_complete_s2) %in% sample_IDs_s2)]
# aqui se conservan solamente los conteos de otus del raton 2 
##---------------------------
meta_samples_phyloseq <- sample_data(meta_samples) # objeto phyloseq

physeq_mice <- merge_phyloseq(otu_table_complete_s2,
                              asv_taxa_simplificada)
physeq_mice2<-phyloseq(physeq_mice, 
                       meta_samples_phyloseq)
