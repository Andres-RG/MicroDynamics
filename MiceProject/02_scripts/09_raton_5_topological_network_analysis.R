# Topological network analysis raton_5 === --- ===
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
##---------------------------
# Se cargan los datos procesados
load(file = "03_out/data/asv_table.RData") # asv_table
load(file = "03_out/data/asv_table_aggregate.RData") # asv_table_cllps
load(file = "03_out/data/basal_period_mouse_5.RData")      # basal_subject5_aggregate
load(file = "03_out/data/fatdiet_period_mouse_5.RData")    # fatdiet_subject5_aggregate
load(file = "03_out/data/recovered1_period_mouse_5.RData") # recover1_subject5_aggregate
load(file = "03_out/data/vancomycin_period_mouse_5.RData") # vancomycin_subject5_aggregate
load(file = "03_out/data/recovered2_period_mouse_5.RData") # recover2_subject5_aggregate
load(file = "03_out/data/gentamicin_period_mouse_5.RData") # gentamicin_subject5_aggregate
load(file = "03_out/data/recovered3_period_mouse_5.RData") # recover3_subject5_aggregate
load(file = "03_out/data/mice_multilayer_network_properties_mouse_5.RData") # mice_ml_properties_s5
##---------------------------
ml_subject5 <- list(
  basal_subject5_aggregate,
  fatdiet_subject5_aggregate,
  recover1_subject5_aggregate,
  vancomycin_subject5_aggregate,
  recover2_subject5_aggregate,
  gentamicin_subject5_aggregate,
  recover3_subject5_aggregate
)
##
# Infiere las co-abundancias a partir de los datos de abundancia
Nbasal_subject5      <- net_inference(basal_subject5_aggregate, "aracne")
Nfatdiet_subject5    <- net_inference(fatdiet_subject5_aggregate, "aracne")
Nrecover1_subject5   <- net_inference(recover1_subject5_aggregate, "aracne")
Nvancomycin_subject5 <- net_inference(vancomycin_subject5_aggregate, "aracne")
Nrecover2_subject5   <- net_inference(recover2_subject5_aggregate, "aracne")
Ngentamicin_subject5 <- net_inference(gentamicin_subject5_aggregate, "aracne")
Nrecover3_subject5   <- net_inference(recover3_subject5_aggregate, "aracne")
##
fatdiet_period_s5    <- list(Nbasal_subject5,
                             Nfatdiet_subject5,
                             Nrecover1_subject5)
vancomyci_period_s5  <- list(Nrecover1_subject5,
                             Nvancomycin_subject5,
                             Nrecover2_subject5)
gentamicin_period_s5 <- list(Nrecover2_subject5,
                             Ngentamicin_subject5,
                             Nrecover3_subject5)
treatments <- c(1,2,3)
##
# Propiedades de las redes multicapa
mice_ml_properties_s5 <- rbind(ml_properties(fatdiet_period_s5, treatments),
                               ml_properties(vancomyci_period_s5, treatments),
                               ml_properties(gentamicin_period_s5, treatments))
##
Stage = c(rep("Fat diet", length(fatdiet_period_s5)),
          rep("Vancomycin", length(vancomyci_period_s5)),
          rep("Gentamicin", length(gentamicin_period_s5))
)
mice_ml_properties_s5 <- cbind(mice_ml_properties_s5, Stage)
## Se genera un objeto que tiene las propiedades de las redes multicapa por 
## cada tratamiento
## Se guarda objeto .RData
# save(mice_ml_properties_s5, file = "03_out/data/mice_multilayer_network_properties_mouse_5.RData")
##---------------------------
names(mice_ml_properties_s5)
# Graficas
## mean degree === === === ==
plot_mean_degree <- ggplot(
  mice_ml_properties_s5,
  aes(x = Treatments,
      y = Mean_degree,
      col = Stage) 
) +
  geom_point() + 
  geom_smooth(
    method = "lm",
    se = F
  ) +
  labs(x = "Period",
       y = "Mean degree",
       tittle = "") +
  theme_minimal() +
  theme(
    # linea del eje
    axis.line = element_line(colour = "black", linewidth = 0.3),
    # eje x
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    # eje y
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    # leyenda
    legend.position = "right",  # Posición de la leyenda
    legend.title = element_text(size = 10, face = "bold"),  # Título de la leyenda
    legend.text = element_text(size = 10),  # Texto de la leyenda
    legend.spacing = unit(0.5, "cm")
  )
## sd degree = === === === ==
plot_sd_degree <- ggplot(
  mice_ml_properties_s5,
  aes(x = Treatments,
      y = sd_degree,
      col = Stage) 
) +
  geom_point() + 
  geom_smooth(
    method = "lm",
    se = F
  ) +
  labs(x = "Period",
       y = "SD of degree",
       tittle = "") +
  theme_minimal() +
  theme(
    # linea del eje
    axis.line = element_line(colour = "black", linewidth = 0.3),
    # eje x
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    # eje y
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    # leyenda
    legend.position = "right",  # Posición de la leyenda
    legend.title = element_text(size = 10, face = "bold"),  # Título de la leyenda
    legend.text = element_text(size = 10),  # Texto de la leyenda
    legend.spacing = unit(0.5, "cm")
  )
## transitivity === === === =
plot_transitivity <- ggplot(
  mice_ml_properties_s5,
  aes(x = Treatments,
      y = Clusterization,
      col = Stage) 
) +
  geom_point() + 
  geom_smooth(
    method = "lm",
    se = F
  ) +
  labs(x = "Period",
       y = "Transitivity",
       tittle = "") +
  theme_minimal() +
  theme(
    # linea del eje
    axis.line = element_line(colour = "black", linewidth = 0.3),
    # eje x
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    # eje y
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    # leyenda
    legend.position = "right",  # Posición de la leyenda
    legend.title = element_text(size = 10, face = "bold"),  # Título de la leyenda
    legend.text = element_text(size = 10),  # Texto de la leyenda
    legend.spacing = unit(0.5, "cm")
  )
## edge density === === === =
plot_edge_density <- ggplot(
  mice_ml_properties_s5,
  aes(x = Treatments,
      y = Edge_density,
      col = Stage) 
) +
  geom_point() + 
  geom_smooth(
    method = "lm",
    se = F
  ) +
  labs(x = "Period",
       y = "Edge density",
       tittle = "") +
  theme_minimal() +
  theme(
    # linea del eje
    axis.line = element_line(colour = "black", linewidth = 0.3),
    # eje x
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    # eje y
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    # leyenda
    legend.position = "right",  # Posición de la leyenda
    legend.title = element_text(size = 10, face = "bold"),  # Título de la leyenda
    legend.text = element_text(size = 10),  # Texto de la leyenda
    legend.spacing = unit(0.5, "cm")
  )
## proportion of linked nodes ===
plot_linked_nodes <- ggplot(
  mice_ml_properties_s5,
  aes(x = Treatments,
      y = Connected_nodes,
      col = Stage) 
) +
  geom_point() + 
  geom_smooth(
    method = "lm",
    se = F
  ) +
  labs(x = "Period",
       y = "Proportion of linked nodes",
       tittle = "") +
  theme_minimal() +
  theme(
    # linea del eje
    axis.line = element_line(colour = "black", linewidth = 0.3),
    # eje x
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    # eje y
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    # leyenda
    legend.position = "right",  # Posición de la leyenda
    legend.title = element_text(size = 10, face = "bold"),  # Título de la leyenda
    legend.text = element_text(size = 10),  # Texto de la leyenda
    legend.spacing = unit(0.5, "cm")
  )
## modularity === === === ===
plot_modularity <- ggplot(
  mice_ml_properties_s5,
  aes(x = Treatments,
      y = Modularity,
      col = Stage) 
) +
  geom_point() + 
  geom_smooth(
    method = "lm",
    se = F
  ) +
  labs(x = "Period",
       y = "Modularity",
       tittle = "") +
  theme_minimal() +
  theme(
    # linea del eje
    axis.line = element_line(colour = "black", linewidth = 0.3),
    # eje x
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    # eje y
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    # leyenda
    legend.position = "right",  # Posición de la leyenda
    legend.title = element_text(size = 10, face = "bold"),  # Título de la leyenda
    legend.text = element_text(size = 10),  # Texto de la leyenda
    legend.spacing = unit(0.5, "cm")
  )
##
# jpeg("03_out/plots/plot_mice_multilayer_network_analysis_mouse_5.jpeg",
#      width = 9333, height = 3700, res = 800, units = "px")
grid.arrange(plot_mean_degree,
             plot_sd_degree,
             plot_transitivity,
             plot_edge_density,
             plot_linked_nodes,
             plot_modularity,
             ncol = 3,
             top = "Multilayer network analysis mouse 5")
# dev.off()