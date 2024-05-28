# For dropdown menu
actionLink <- function(inputId, ...) {
  tags$a(
    href = "javascript:void",
    id = inputId,
    class = "action-button",
    ...
  )
}


### **Function for creating networks that the distance between nodes is 1 HD**

# Create a net link function using 1HD as distance
simple_net <- function(x) {
  uniqueids2 <- unique(as.character(x$junction))
  links <- stringdistmatrix(uniqueids2, uniqueids2, method = "hamming")
  rownames(links) <- uniqueids2
  colnames(links) <- uniqueids2
  net <- graph_from_adjacency_matrix(links < 2, diag = F, mode = "undirected")
}


### Function for generating key network parameters for clusters

generate_network_parameters <- function(input_data) {
  # Remove singleton clusters (only 1 sequence inside)
  # data_pars <- input_data %>%
  #   group_by(cluster) %>%
  #   dplyr::mutate(samples = n_distinct(sample)) %>%
  #   dplyr::mutate(seqs = n_distinct(junction)) %>%
  #   filter(seqs > 1)
  
  data_pars <- input_data %>%
    group_by(cluster) %>%
    dplyr::mutate(seqs = n_distinct(junction)) %>%
    filter(seqs > 1)
  
  
  # Create the pivot table of clonotypes and samples
  # data <- data_pars %>%
  #   distinct(cluster, junction, sample) %>%
  #   group_by(cluster, sample, junction) %>%
  #   dplyr::summarise(count = n(), .groups = "drop") %>%
  #   pivot_wider(names_from = sample, values_from = count, values_fill = 0)
  
  data <- data_pars %>% distinct(cluster, junction)
  
  # Get the actual connection nodes with 1 HD
  pw_edge <- data %>%
    group_by(cluster) %>%
    do(as.data.frame(get.edgelist(simple_net(.)))) %>%
    as_tibble()
  
  # Get the maximum HD in the clusters
  pw_edge2 <- data %>%
    group_by(cluster) %>%
    reframe(hds = unlist(as.list(stringdistmatrix(unique(as.character(junction)), unique(as.character(junction)), method = "hamming")))) %>%
    group_by(cluster) %>%
    mutate("max_HD" = max(hds)) %>%
    select(-hds) %>%
    distinct()
  
  # Construct the net using unique junction sequences for plot
  list_of_edges <- split(pw_edge, pw_edge$cluster)
  
  # Creating a list of igraph objects based on each subset of edges
  net <- lapply(list_of_edges, function(edges) {
    graph_from_data_frame(d = edges[, 2:3], vertices = unique(c(edges$V1, edges$V2)), directed = FALSE)
  })
  
  
  # Get the network parameters from igraph function
  diameters <- sapply(net, igraph::diameter)
  densitys <- sapply(net, igraph::edge_density)
  mean_distances <- sapply(net, igraph::mean_distance)
  ecount <- sapply(net, igraph::ecount)
  
  degree <- sapply(net, function(g) max(igraph::degree(g)))
  mean_degrees <- sapply(net, function(g) mean(igraph::degree(g)))
  
  degree_data <- do.call(rbind, lapply(1:length(net), function(i) {
    data.frame(cluster = names(net[i]), n = igraph::degree(net[[i]]))
  }))
  
  clusters <- names(list_of_edges)
  
  # Create a dataframe with cluster, diameter, and mean distance
  result_df <- data.frame(
    cluster = clusters,
    diameter = diameters,
    density = densitys,
    mean_distance = mean_distances,
    ecount = ecount,
    degree = degree,
    mean_degree = mean_degrees,
    spread = ecount / (degree * mean_degrees)
  )
  
  result_df$cluster <- as.numeric(result_df$cluster)
  # Add max HD to the result
  result_df <- inner_join(result_df, pw_edge2, by = "cluster")
  
  # Collect all the parameters and create the final dataframe
  data_pars <- inner_join(result_df, data_pars, by = "cluster") %>%
    mutate(cluster_density = seqs / max_HD) %>%
    mutate_if(is.numeric, round, digits = 2)
}


############################################

### Function for creating network layout (fr, kk, kk3)

############################################


# Function to plot a network visualization based on the provided datasets and parameters.
# Parameters:
# - test: A small dataset with cluster and locus information.
# - data: A full dataset with clonotypes, locus, cluster, color settings, compartment or diagnosis.
# - group: Determines which feature the node should represent in the network.
# - color_col: Specifies the colors for the group presentation.
# - scale: Factor to scale the size of nodes in the network plot.
# - circle_color: Color for the frame of the network layout.
# - layout: Specifies the layout algorithm for the network visualization. Options include
#   'fr' (Fruchterman-Reingold), 'kk2' (Kamada-Kawai 2D), 'kk3' (Kamada-Kawai 3D),
#   and 'drl' (drL for large-scale graphs).
# - legend: 'on' to display the group information as a legend in the plot.
# - seed: A random seed number for adjusting the layout of network plot

plot_network <- function(test, data, group, color_col, scale, circle_color, layout, legend, seed, text_size) {
  # Transform the group and color column names into symbols for dynamic referencing.
  
  group_sym <- sym(group)
  color_sym <- sym(color_col)
  
  # Pivot the data for network analysis and create a semi-join with the 'test' dataset.
  # This process selects distinct clusters and junctions, groups them, and prepares them for plotting.
  data_pivot <- data %>%
    semi_join(test, by = c("cluster", "locus")) %>%
    distinct(cluster, junction, !!group_sym) %>%
    group_by(cluster, !!group_sym, junction) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = !!group_sym, values_from = count, values_fill = 0) %>%
    arrange(cluster, junction)
  
  
  # Create a pairwise edge list for network analysis and visualization.
  pw_edge <- data_pivot %>%
    group_by(cluster) %>%
    do(as.data.frame(get.edgelist(simple_net(.)))) %>%
    as_tibble()
  
  # Construct the net using unique junction sequences for plot
  list_of_edges <- split(pw_edge, pw_edge$cluster)
  
  # Create a list of igraph objects based on each subset of edges
  net <- lapply(list_of_edges, function(edges) {
    graph_from_data_frame(d = edges[, 2:3], vertices = unique(c(edges$V1, edges$V2)), directed = FALSE)
  })
  
  
  # Iteratively bind data frames to prepare them for visualization.
  out <- NULL
  for (i in seq_along(net)) {
    for (vrt in 1:length(V(net[[i]]))) {
      vertex_in_net <- data_pivot[which(data_pivot$junction == V(net[[i]])$name[vrt])[1], ]
      out <- rbind(out, vertex_in_net)
    }
  }
  
  out2 <- NULL
  for (i in seq_along(net)) {
    for (vrt in 1:length(V(net[[i]]))) {
      # vertex_in_net <- data[which(data$junction==V(net[[i]])$name[vrt]),]
      vertex_in_net <- data[which(data$junction == V(net[[i]])$name[vrt])[1], ]
      out2 <- rbind(out2, vertex_in_net)
    }
  }
  
  
  
  if (!"reads" %in% names(out2)) {
    # If 'reads' does not exist, create it and set all values to 1
    out2$reads <- 1000
  } else if  (any(is.na(out2$reads) | out2$reads == "NA")) {
    out2$reads <- 1000
  } 
  
  
  
  # Group data and calculate the total reads for each junction and cluster.
  out2 <- out2 %>% group_by(cluster, junction, !!group_sym) %>%
    mutate(all_reads=sum(as.numeric(reads)))
  
  
  numeric_data <- (log1p(as.numeric(out2$all_reads)))
  
  # Prepare data for determining the size of pie charts in the network nodes.
  am <- out %>%
    select((match("junction", names(out)) + 1):ncol(out)) %>%
    simplify2array(.)
  
  # Note the pie nodes were determined by the distinct number of
  i <- 1
  value <- lapply(seq_len(nrow(am)), function(i) as.numeric(unlist(am[i, ])))
  value <- lapply(value, function(x) {
    x[x != 0]
  })
  
  # Modify the network object names and prepare for merging.
  nets_modified <- lapply(seq_along(net), function(i) {
    name <- names(net)[i]
    G <- net[[i]]
    V(G)$name <- paste(V(G)$name, name, sep = "_")
    return(G)
  })
  
  # Merge all individual networks into a single disjoint network.
  merged_net <- do.call(disjoint_union, nets_modified)
  
  # Set layout parameters based on the user's choice of layout algorithm.
  set.seed(seed)
  if (layout == "fr") {
    lo <- layout_with_fr(merged_net, dim = 2, grid = "nogrid", niter = 1500)
  } else if (layout == "kk2") {
    lo <- layout_with_kk(merged_net, dim = 2)
  } else if (layout == "kk3") {
    lo <- layout_with_kk(merged_net, dim = 3)
  } else {
    lo <- layout_with_drl(merged_net, dim = 2)
  }
  
  # Normalize coordinates and shift layout to center at (0,0).
  lo <- norm_coords(lo)
  
  # Calculate mean x and y coordinates
  mean_x <- mean(lo[, 1])
  mean_y <- mean(lo[, 2])
  
  # Shift layout to be centered at (0,0)
  lo[, 1] <- lo[, 1] - mean_x
  lo[, 2] <- lo[, 2] - mean_y
  
  
  max_dist <- max(apply(lo, 1, function(coord) sqrt(sum(coord^2))))
  adjusted_max_dist <- 1.1 * max_dist
  
  # Calculate the relative size of the nodes based on max_dist
  relative_size <- numeric_data / max(numeric_data) * adjusted_max_dist
  
  # Scale these sizes to fit the plotting region
  scaled_size <- (relative_size / adjusted_max_dist) * par("din")[1] * scale
  
  
  # Group by junction and count the number of colors for each junction
  junction_groups <- data %>%
    semi_join(test, by = c("cluster", "locus")) %>%
    group_by(cluster, junction) %>%
    dplyr::mutate(color_count = n_distinct(!!color_sym)) %>%
    mutate(jc = paste(junction, cluster, sep = "_")) %>%
    distinct(cluster, junction, color_count, jc)
  
  # Identify pie nodes and normal nodes
  pie_nodes <- (junction_groups$jc[junction_groups$color_count > 1])
  normal_nodes <- (junction_groups$jc[junction_groups$color_count == 1])
  
  
  # Create the same format of column (junction_cluster) to match the network node format
  sub <- data %>%
    filter(cluster %in% test$cluster) %>%
    distinct(cluster, junction, !!color_sym, !!group_sym)
  sub$jc <- paste(sub$junction, sub$cluster, sep = "_")
  
  # Correspond the color and the node
  split_colors <- split(sub[[color_col]], sub$jc)
  
  new_split_colors <- lapply(V(merged_net)$name, function(node) {
    if (node %in% names(split_colors)) {
      return(split_colors[[node]])
    } else {
      return(NA)
    }
  })
  
  # Match the node and colors
  names(new_split_colors) <- V(merged_net)$name
  split_colors <- lapply(new_split_colors, unique)
  
  # Prepare colors for pie and regular nodes based on the network data.
  pie_colors <- vector("list", length = vcount(merged_net))
  regular_colors <- rep(NA, vcount(merged_net))
  
  # Add colors to the node based on pie or regular
  for (node in V(merged_net)$name) {
    idx <- which(V(merged_net)$name == node)
    if (length(idx) > 0) {
      idx <- idx[1] # Ensure a single index is used
      if (node %in% pie_nodes) {
        if (!is.null(split_colors[[node]])) {
          pie_colors[[idx]] <- split_colors[[node]]
        }
      } else if (node %in% normal_nodes) {
        if (!is.null(split_colors[[node]])) {
          regular_colors[idx] <- split_colors[[node]][1]
        }
      }
    }
  }
  
  
  # Assign shapes and frame colors to the vertices in the network.
  vertex_shapes <- rep("circle", vcount(merged_net))
  vertex_shapes[which(V(merged_net)$name %in% pie_nodes)] <- "pie"
  
  vertex_frame_colors <- rep("black", vcount(merged_net))
  
  # Plot the network with the specified parameters, including colors, shapes, and layout.
  plot(merged_net,
       rescale = FALSE,
       axes = FALSE,
       asp = 1,
       xlim = c(-adjusted_max_dist, adjusted_max_dist),
       ylim = c(-adjusted_max_dist, adjusted_max_dist),
       vertex.color = regular_colors,
       vertex.pie.color = pie_colors,
       vertex.frame.color = vertex_frame_colors,
       vertex.shape = vertex_shapes,
       layout = lo,
       edge.label = NA,
       vertex.pie = value,
       edge.color = "darkgrey",
       vertex.pie.lty = "blank",
       vertex.pie.width = 0.35,
       vertex.frame.width = 0.2,
       edge.width = 0.6, vertex.label = NA,
       vertex.size = scaled_size
  )
  
  # Add a circle using the circle function
  draw.circle(0, 0, adjusted_max_dist, border = circle_color, lwd = 1)
  
  # If the legend is enabled, create and display a legend based on unique group-color combinations.
  if (legend == "on") {
    # Get unique combinations
    unique_combinations <- unique(sub[, c(group, color_col)])
    
    # Set names based on levels
    v_colors <- setNames(unique_combinations[[color_col]], unique_combinations[[group]])
    
    legend("bottom",
           bty = "n", legend = names(v_colors),
           fill = v_colors, ncol = 2, inset = c(0, -.005),
           y.intersp = 0.7, cex = 0.9, text.font = text_size
    )
  }
}


############################################

### Function for interactive network layout I

############################################

plot_network_vis1<-function(test, data){
  
  # Transform the group and color column names into symbols for dynamic referencing.
  
  
  # Pivot the data for network analysis and create a semi-join with the 'test' dataset.
  # This process selects distinct clusters and junctions, groups them, and prepares them for plotting.
  data_pivot <- data %>% 
    semi_join(test, by = c("cluster", "locus")) %>%
    distinct(cluster, junction) 
  
  
  # Create a pairwise edge list for network analysis and visualization.
  pw_edge<-data_pivot %>% group_by(cluster) %>% 
    do(as.data.frame(get.edgelist(simple_net(.)))) %>%
    as_tibble() 
  
  # Construct the net using unique junction sequences for plot
  list_of_edges <- split(pw_edge, pw_edge$cluster)
  
  # Create a list of igraph objects based on each subset of edges
  net <- lapply(list_of_edges, function(edges) {
    graph_from_data_frame(d=edges[,2:3], vertices=unique(c(edges$V1, edges$V2)), directed=FALSE)
  })
  
  
  # Modify the network object names and prepare for merging.   
  nets_modified <- lapply(seq_along(net), function(i) {
    name <- names(net)[i]
    G<-net[[i]]
    V(G)$name <- paste(V(G)$name, name, sep = "_")
    return(G)
  })
  
  # Merge all individual networks into a single disjoint network. 
  merged_net <- do.call(disjoint_union, nets_modified)  
  
  
  
  # Extract nodes
  nodes_name <-data.frame(
    id = as.numeric(V(merged_net)), full_names = V(merged_net)$name) %>% 
    separate(full_names, into = c("label", "cluster"), sep = "_") %>%
    select(id, label)
  
  rownames(nodes_name) <- NULL
  
  # Extract edges
  edges <- data.frame(
    from = gsub("_.*[0-9]","",ends(merged_net, E(merged_net)) [, 1]),
    to = gsub("_.*[0-9]","",ends(merged_net, E(merged_net)) [, 2])
    # ... include other edge attributes if necessary
  )
  
  # Create a lookup table for labels to ids
  lookup_table <- nodes_name[, c("label", "id")]
  
  
  # Replace 'from' labels with ids
  edges <- merge(edges, lookup_table, by.x = "from", by.y = "label")
  
  edges<-edges %>% select(-from) %>% dplyr::rename(from = id)
  
  
  # Replace 'to' labels with ids
  edges <- merge(edges, lookup_table, by.x = "to", by.y = "label")
  edges<- edges %>% select(-to) %>% dplyr::rename(to = id)
  edges$from <- as.numeric(edges$from)
  edges$to <- as.numeric(edges$to)
  
  
  visNetwork(nodes = nodes_name,edges) %>% 
    visOptions(highlightNearest=TRUE, 
               nodesIdSelection = TRUE)  %>%      #allow for long click to select additional nodes
    visInteraction(multiselect = TRUE)
  
  
}


############################################

### Function for interactive network layout II

############################################


plot_network_vis2<-function(test, data, layout, group, color_col){
  
  
  # Transform the group and color column names into symbols for dynamic referencing.
  
  group_sym <- sym(group)
  color_sym <- sym(color_col)
  
  # Pivot the data for network analysis and create a semi-join with the 'test' dataset.
  # This process selects distinct clusters and junctions, groups them, and prepares them for plotting.
  data_pivot <- data %>%
    semi_join(test, by = c("cluster", "locus")) %>%
    distinct(cluster, junction, !!group_sym) %>%
    group_by(cluster, !!group_sym, junction) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = !!group_sym, values_from = count, values_fill = 0) %>%
    arrange(cluster, junction)
  
  
  # Create a pairwise edge list for network analysis and visualization.
  pw_edge <- data_pivot %>%
    group_by(cluster) %>%
    do(as.data.frame(get.edgelist(simple_net(.)))) %>%
    as_tibble()
  
  # Construct the net using unique junction sequences for plot
  list_of_edges <- split(pw_edge, pw_edge$cluster)
  
  # Create a list of igraph objects based on each subset of edges
  net <- lapply(list_of_edges, function(edges) {
    graph_from_data_frame(d = edges[, 2:3], vertices = unique(c(edges$V1, edges$V2)), directed = FALSE)
  })
  
  
  # Iteratively bind data frames to prepare them for visualization.
  out <- NULL
  for (i in seq_along(net)) {
    for (vrt in 1:length(V(net[[i]]))) {
      vertex_in_net <- data_pivot[which(data_pivot$junction == V(net[[i]])$name[vrt])[1], ]
      out <- rbind(out, vertex_in_net)
    }
  }
  
  out2 <- NULL
  for (i in seq_along(net)) {
    for (vrt in 1:length(V(net[[i]]))) {
      # vertex_in_net <- data[which(data$junction==V(net[[i]])$name[vrt]),]
      vertex_in_net <- data[which(data$junction == V(net[[i]])$name[vrt])[1], ]
      out2 <- rbind(out2, vertex_in_net)
    }
  }
  
  
  
  if (!"reads" %in% names(out2)) {
    # If 'reads' does not exist, create it and set all values to 1
    out2$reads <- 1000
  } else if  (any(is.na(out2$reads) | out2$reads == "NA")) {
    out2$reads <- 1000
  } 
  
  
  
  # Group data and calculate the total reads for each junction and cluster.
  out2 <- out2 %>% group_by(cluster, junction, !!group_sym) %>% 
    mutate(all_reads=sum(as.numeric(reads)))
  
  
  numeric_data <- (log1p(as.numeric(out2$all_reads)))
  
  # Prepare data for determining the size of pie charts in the network nodes.
  am <- out %>%
    select((match("junction", names(out)) + 1):ncol(out)) %>%
    simplify2array(.)
  
  # Note the pie nodes were determined by the distinct number of
  i <- 1
  value <- lapply(seq_len(nrow(am)), function(i) as.numeric(unlist(am[i, ])))
  value <- lapply(value, function(x) {
    x[x != 0]
  })
  
  # Modify the network object names and prepare for merging.
  nets_modified <- lapply(seq_along(net), function(i) {
    name <- names(net)[i]
    G <- net[[i]]
    V(G)$name <- paste(V(G)$name, name, sep = "_")
    return(G)
  })
  
  # Merge all individual networks into a single disjoint network.
  merged_net <- do.call(disjoint_union, nets_modified)
  
  
  # Extract nodes
  nodes_name <-data.frame(
    id = as.numeric(V(merged_net)), full_names = V(merged_net)$name) %>% 
    separate(full_names, into = c("label", "cluster"), sep = "_") %>%
    select(id, label)
  
  rownames(nodes_name) <- NULL
  
  # Extract edges
  edges <- data.frame(
    from = gsub("_.*[0-9]","",ends(merged_net, E(merged_net)) [, 1]),
    to = gsub("_.*[0-9]","",ends(merged_net, E(merged_net)) [, 2])
  )
  
  
  karate <- igraph::graph_from_data_frame(edges, 
                                          directed = FALSE)
  
  for(node in V(karate)$name) {
    # Find the matching row in the out2 dataframe
    matching_row <- which(out2$junction == node)
    
    # Check if a matching row is found
    if(length(matching_row) > 0) {
      # Get the 'cat' value from the matching row
      group_values <- out2[[color_sym]][matching_row]
      unique_group_values <- unique(sort(group_values))
      # Assign the 'cat' value to the 'group' attribute of the node
      if (length(unique_group_values) == 1) {
        V(karate)$color[V(karate)$name == node] <- unique_group_values
        V(karate)$shape[V(karate)$name == node] <- "circle"
      } else {
        V(karate)$color[V(karate)$name == node] <- "black"
        V(karate)$shape[V(karate)$name == node] <- "square" 
      }
      
    }
  }
  
  visNetwork::visIgraph(karate, layout = layout, 
                        randomSeed = 123) %>% visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
  
}

############################################

### Function for plotting seqlogo for each cluster

############################################


plot_seqlogo <- function(test, data_ori, title_size) {
  plot_list <- list()
  for (j in 1:length(as.list(unique(test$cluster)))) {
    # Prepare the data with v, j gene subgroup information and actural junction sequences
    data <- data_ori %>%
      semi_join(test[test$cluster == test$cluster[j], ], by = c("locus", "cluster")) %>%
      distinct(cluster, junction, v_call, j_call) %>%
      mutate(v_subgroup = gsub("-.*[0-9]", "", v_call))
    
    if (length(data$junction) > 2) {
      # Sequences alignment
      plot_list[[j]] <- ggplot() +
        geom_logo(data = as.character(AlignSeqs(AAStringSet(data$junction)))) +
        ggtitle(paste0("Cluster ID: ", data$cluster, ", Cluster size: ", (unique(test$seqs)), ", V subgroup: ", data$v_subgroup, ", J subgroup: ", data$j_call)) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white")
        ) +
        theme(
          axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size=title_size ),
        ) +
        theme(legend.position = "none")
    } else {
      print("Seqlogo function needs at least 3  sequences ")
    }
  }
  
  gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
}