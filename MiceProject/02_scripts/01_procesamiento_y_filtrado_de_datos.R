# Librerias necesarias
library(readr)
library(mlBioNets)
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/asv_table.RData") # asv_table
load(file = "03_out/data/asv_table_aggregate.RData") # asv_table_cllps
load(file = "03_out/data/basal_period_mouse_2.RData") # basal_subject2_aggregate
load(file = "03_out/data/fatdiet_period_mouse_2.RData") # fatdiet_subject2_aggregate
load(file = "03_out/data/recovered1_period_mouse_2.RData") # recover1_subject2_aggregate
load(file = "03_out/data/vancomycin_period_mouse_2.RData") # vancomycin_subject2_aggregate
load(file = "03_out/data/recovered2_period_mouse_2.RData") # recover2_subject2_aggregate
load(file = "03_out/data/gentamicin_period_mouse_2.RData") # gentamicin_subject2_aggregate
load(file = "03_out/data/recovered3_period_mouse_2.RData") # recover3_subject2_aggregate
##---------------------------
# Se cargan los datos en el entorno de R
# asv abundance table
asv_table <- as.data.frame(read_tsv("01_raw_data/counts.tsv"))
asv_ids <- as.vector(asv_table[,1])
asv_table <- asv_table[,-1]
rownames(asv_table) <- asv_ids
#   Se puso el nombre de los asv en los renglones de la tabla, ya que venían 
#   agregados como una columna extra
#   Se guarda como archivo .RData
# save(asv_table, file = "03_out/data/asv_table.RData")
# asv metadata table
asv_meta <- as.data.frame(read_tsv("01_raw_data/metadata.tsv"))
# perturbations table
asv_pert <- as.data.frame(read_tsv("01_raw_data/perturbations.tsv"))
# taxonomic table
asv_taxa<-as.data.frame(read_tsv("01_raw_data/rdp_species.tsv"))
# aggregate table at genus level
asv_table_aggregate <- T_collapse(is_phyloseq = F,
                                  T_table = asv_taxa,
                                  O_table = asv_table,
                                  names_level = "Genus")
#   Junta todas las abundancias de los taxa y los guarda en un objeto, 
#   con todas las series de tiempo
#   Se guarda el objeto como .RData
# save(asv_table_aggregate, file = "03_out/data/asv_table_aggregate.RData")
##---------------------------
# Se separa la tabla por cada subject, que es cada raton
# el raton 1 no tenpia perturbaciones, se toman unicamente los ratones 
# 2, 3, 4 y 5
## RATON 2
subject2 <- asv_meta[which(asv_meta$subject == 2),]
subject2 <- subject2[order(subject2$time), ]
## RATON 3
subject3 <- asv_meta[which(asv_meta$subject == 3),]
subject3 <- subject3[order(subject3$time), ]
## RATON 4
subject4 <- asv_meta[which(asv_meta$subject == 4),]
subject4 <- subject4[order(subject4$time), ]
## RATON 5
subject5 <- asv_meta[which(asv_meta$subject == 5),]
subject5 <- subject5[order(subject5$time), ]
##---------------------------