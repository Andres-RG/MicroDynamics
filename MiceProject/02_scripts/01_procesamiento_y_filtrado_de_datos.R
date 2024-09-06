# Librerias necesarias
library(readr)

##-------


# asv abundance table
asv_table <- as.data.frame(read_tsv("01_raw_data/counts.tsv"))
asv_ids <- as.vector(asv_table[,1])
asv_table <- asv_table[,-1]
rownames(asv_table) <- asv_ids
# Se puso el nombre de los asv en los renglones de la tabla, ya que venían 
# agregados como una columna extra
# Se guarda como archivo .RData
save(asv_table, file = "03_out/data/asv_table.RData")

# asv metadata table
asv_meta <- as.data.frame(read_tsv("01_raw_data/metadata.tsv"))

# perturbations table
asv_pert <- as.data.frame(read_tsv("01_raw_data/perturbations.tsv"))

# taxonomic table
