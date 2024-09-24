# Topological network analysis raton_2 === --- ===
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
load(file = "03_out/data/basal_period_mouse_2.RData")      # basal_subject2_aggregate
load(file = "03_out/data/fatdiet_period_mouse_2.RData")    # fatdiet_subject2_aggregate
load(file = "03_out/data/recovered1_period_mouse_2.RData") # recover1_subject2_aggregate
load(file = "03_out/data/vancomycin_period_mouse_2.RData") # vancomycin_subject2_aggregate
load(file = "03_out/data/recovered2_period_mouse_2.RData") # recover2_subject2_aggregate
load(file = "03_out/data/gentamicin_period_mouse_2.RData") # gentamicin_subject2_aggregate
load(file = "03_out/data/recovered3_period_mouse_2.RData") # recover3_subject2_aggregate
load(file = "03_out/data/mice_multilayer_network_properties_mouse_2.RData") # mice_ml_properties
##---------------------------
ml_subject2 <- list(
  basal_subject2_aggregate,
  fatdiet_subject2_aggregate,
  recover1_subject2_aggregate,
  vancomycin_subject2_aggregate,
  recover2_subject2_aggregate,
  gentamicin_subject2_aggregate,
  recover3_subject2_aggregate
)
##
# Infiere las co-abundancias a partir de los datos de abundancia
Nbasal_subject2      <- net_inference(basal_subject2_aggregate, "aracne")
Nfatdiet_subject2    <- net_inference(fatdiet_subject2_aggregate, "aracne")
Nrecover1_subject2   <- net_inference(recover1_subject2_aggregate, "aracne")
Nvancomycin_subject2 <- net_inference(vancomycin_subject2_aggregate, "aracne")
Nrecover2_subject2   <- net_inference(recover2_subject2_aggregate, "aracne")
Ngentamicin_subject2 <- net_inference(gentamicin_subject2_aggregate, "aracne")
Nrecover3_subject2   <- net_inference(recover3_subject2_aggregate, "aracne")
##
fatdiet_period_s2    <- list(Nbasal_subject2,
                             Nfatdiet_subject2,
                             Nrecover1_subject2)
vancomyci_period_s2  <- list(Nrecover1_subject2,
                             Nvancomycin_subject2,
                             Nrecover2_subject2)
gentamicin_period_s2 <- list(Nrecover2_subject2,
                             Ngentamicin_subject2,
                             Nrecover3_subject2)
treatments <- c(1,2,3)
##
# Propiedades de las redes multicapa
mice_ml_properties <- rbind(ml_properties(fatdiet_period_s2, treatments),
                            ml_properties(vancomyci_period_s2, treatments),
                            ml_properties(gentamicin_period_s2, treatments))
##
Stage = c(rep("Fat diet", length(fatdiet_period_s2)),
          rep("Vancomycin", length(vancomyci_period_s2)),
          rep("Gentamicin", length(gentamicin_period_s2))
          )
mice_ml_properties <- cbind(mice_ml_properties, Stage)
## Se genera un objeto que tiene las propiedades de las redes multicapa por 
## cada tratamiento
## Se guarda objeto .RData
# save(mice_ml_properties, file = "03_out/data/mice_multilayer_network_properties_mouse_2.RData")
##---------------------------
names(mice_ml_properties)
# Graficas
## mean degree === === === ==
plot_mean_degree <- ggplot(
  mice_ml_properties,
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
  mice_ml_properties,
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
  mice_ml_properties,
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
  mice_ml_properties,
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
  mice_ml_properties,
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
  mice_ml_properties,
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
# jpeg("03_out/plots/plot_mice_multilayer_network_analysis_mouse_2.jpeg",
#      width = 9333, height = 3700, res = 800, units = "px")
grid.arrange(plot_mean_degree,
             plot_sd_degree,
             plot_transitivity,
             plot_edge_density,
             plot_linked_nodes,
             plot_modularity,
             ncol = 3,
             top = "Multilayer network analysis")
# dev.off()