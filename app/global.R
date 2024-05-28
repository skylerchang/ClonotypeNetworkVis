library(shiny)
library(vroom)
library(grid)
library(dplyr)
library(ggplot2)
library(reshape)
library(forcats)
library(ggbeeswarm)
library(tidyr)
library(stringr)
library(visNetwork)
library(igraph)
library(stringdist)
library(RColorBrewer)
library(plotrix)
library(ggrepel)
library(beeswarm)
library(DECIPHER)
library(ggseqlogo)
library(Biostrings)
library(ggpubr)
library(gridExtra)
library(UpSetR)
library(shinycssloaders)
library(shinythemes)
library(ggvis)
library(plotly)
library(shinydashboard)
library(reticulate)
library(DT)
library(processx)
library(htmlwidgets)
library(webshot2)
library(ape)
library(phangorn)
library(treeio)
library(msa)
library(ggtree)
library(here)




# increase the uploading file size limit to 2000M
options(shiny.maxRequestSize = 2000 * 1024^2)

button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;

/* Change the text size to 15 pixels. */
font-size: 15px;
}"


# Color selection
color_vars <- c(
  "Set1" = "Set1",
  "Accent" = "Accent",
  "Paired" = "Paired",
  "Dark2" = "Dark2",
  "Set2" = "Set2",
  "Set3" = "Set3"
)


axis_vars <- c(
  "Diameter" = "diameter",
  "Density" = "density",
  "Mean distance" = "mean_distance",
  "Edge number" = "ecount",
  "Cluster centrality" = "degree",
  "Mean centrality" = "mean_degree",
  "Spread level" = "spread",
  "Max Hamming distance (HD)" = "max_HD",
  "Junctional length" = "junction_length",
  "Clonotype number" = "seqs",
  "Cluster density" = "cluster_density"
)



layout_vars <- c(
  "layout_with_fr" = "layout_with_fr",
  "layout_with_drl" = "layout_with_drl",
  "layout_with_kk" = "layout_with_kk",
  "layout_with_dh" = "layout_with_dh",
  "layout_with_gem" = "layout_with_gem",
  "layout_with_graphopt" = "layout_with_graphopt",
  "layout_with_lgl" = "layout_with_lgl",
  "layout_with_mds" = "layout_with_mds",
  "layout_with_sugiyama" = "layout_with_mds"
)

layouts <- c("fr" = "fr", "kk2" = "kk2", "kk3" = "kk3")

ggtree_layouts <- c(
  "roundrect" = "roundrect",
  "slanted" = "slanted",
  "ellipse" = "ellipse",
  "circular" = "circular",
  "equal_angle" = "equal_angle",
  "daylight" = "daylight"
)

source(here("R", "functions.R"))


