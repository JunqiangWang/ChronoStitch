


#'
#'
#'
CalculateMST <- function(seu,
                         reduction = "umap",
                         dist.method = "euclidean") {
  # Load required libraries
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  library(Matrix)
  library(igraph)
  library(dplyr)
  
  # Extract embedding from the specified dimensionality reduction
  embedding <- Embeddings(seu, reduction = reduction)
  
  # Compute distance matrix
  dist_obj <- dist(embedding, method = dist.method)
  dist_mat <- as.matrix(dist_obj)
  
  # Assign consistent row and column names
  cell_names <- colnames(seu)
  rownames(dist_mat) <- cell_names
  colnames(dist_mat) <- cell_names
  
  # Convert to sparse matrix if not already
  if (!inherits(dist_mat, "dgCMatrix")) {
    dist_mat <- Matrix(dist_mat, sparse = TRUE)
  }
  
  # Convert sparse matrix to edge list (triplet format)
  triplet_df <- as.data.frame(summary(dist_mat))
  colnames(triplet_df) <- c("from", "to", "weight")
  
  # Remove self-loops and duplicate edges
  edge_df <- triplet_df %>%
    filter(from != to) %>%
    mutate(edge_id = paste0(pmin(from, to), "_", pmax(from, to))) %>%
    distinct(edge_id, .keep_all = TRUE) %>%
    select(from, to, weight)
  
  # Create undirected graph
  g <- graph_from_data_frame(edge_df, directed = FALSE)
  
  # Compute Minimum Spanning Tree (MST)
  mst_graph <- igraph::mst(g, weights = E(g)$weight)
  
  V(mst_graph)$name <- colnames(seu)
  
  seu@misc$mst_graph<-mst_graph
  
  return(seu)
  
}




#' 
#' 
#' 
PlotMST<-function(seu,
                  reduction="umap" ){
  
  mst_graph<-seu@misc$mst_graph
  
  V(mst_graph)$name <- colnames(seu) # from earlier
  
  # ---- 10. Plot MST ----
  # plot(mst_graph,
  #      vertex.label = NA,
  #      vertex.size = 5,
  #      edge.width = 1,
  #      main = "MST")
  
  # ---- 10. Plot MST over UMAP ----
  # Get 2D layout from Seurat (UMAP)
  umap_layout <- Embeddings(seu, reduction = reduction)
  
  # Match node names in graph to UMAP rownames (which should be cell names)
  # Ensure correct node labeling
  # Extract layout in correct order (igraph uses V(mst_knn)$name)
  layout_coords <- umap_layout[V(mst_graph)$name, ]
  
  # Plot MST using UMAP layout as coordinates
  plot(mst_graph,
       layout = layout_coords,
       vertex.label = NA,
       vertex.size = 2,
       edge.width = 1,
       edge.color = "grey30",
       vertex.color = "steelblue",
       main = "MST over UMAP")
  
}



#'
#'
#' 
GetRootCell <- function(seu,
                        group_var = NULL,
                        root_cell_type = NULL,
                        return_root_cell = FALSE
) {
  # Load required libraries
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' is required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  
  library(igraph)
  library(FNN)
  library(Matrix)
  
  # ---- 1. Validate inputs ----
  if (is.null(group_var) || is.null(root_cell_type)) {
    stop("Both 'group_var' and 'root_cell_type' must be specified.")
  }
  
  if (!group_var %in% colnames(seu@meta.data)) {
    stop(paste0("'", group_var, "' not found in Seurat metadata."))
  }
  
  if (is.null(seu@misc$mst)) {
    stop("No MST found in seu@misc. Run CalculateMST() first.")
  }
  
  mst <- seu@misc$mst
  
  # ---- 2. Get UMAP coordinates ----
  umap_coords <- Embeddings(seu, "umap")
  umap_coords <- umap_coords[colnames(seu), , drop = FALSE]
  
  # ---- 3. Subset cells of target type ----
  cell_ids <- rownames(seu@meta.data[seu@meta.data[[group_var]] %in% root_cell_type, ])
  if (length(cell_ids) == 0) {
    stop("No cells found for the specified root_cell_type.")
  }
  
  cell_coords <- umap_coords[cell_ids, , drop = FALSE]
  
  # ---- 4. Extract MST vertex names and layout ----
  mst_vertex_names <- V(mst)$name
  mst_layout <- umap_coords[mst_vertex_names, , drop = FALSE]
  
  # ---- 5. Map each cell to its nearest MST vertex using FNN ----
  knn_map <- get.knnx(mst_layout, cell_coords, k = 1)
  closest_vertex_indices <- knn_map$nn.index[, 1]
  closest_vertex_names <- mst_vertex_names[closest_vertex_indices]
  
  # ---- 6. Find the most common vertex as root ----
  root_vertex <- names(which.max(table(closest_vertex_names)))
  cat("Most common root MST vertex for", root_cell_type, ":", root_vertex, "\n")
  
  
  if(isTRUE(return_root_cell)){
    return(root_vertex)
  }else{
    seu@misc$root_cell<-root_vertex
  }
  
}




GetRootCell <- function(seu,
                        group_var = NULL,
                        root_cell_type = NULL,
                        return_root_cell = FALSE
) {
  # Load required libraries
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' is required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  
  library(igraph)
  library(FNN)
  library(Matrix)
  
  # ---- 1. Validate inputs ----
  if (is.null(group_var) || is.null(root_cell_type)) {
    stop("Both 'group_var' and 'root_cell_type' must be specified.")
  }
  
  if (!group_var %in% colnames(seu@meta.data)) {
    stop(paste0("'", group_var, "' not found in Seurat metadata."))
  }
  
  if (is.null(seu@misc$mst)) {
    stop("No MST found in seu@misc. Run CalculateMST() first.")
  }
  
  mst <- seu@misc$mst
  
  # ---- 2. Get UMAP coordinates ----
  umap_coords <- Embeddings(seu, "umap")
  umap_coords <- umap_coords[colnames(seu), , drop = FALSE]
  
  # ---- 3. Subset cells of target type ----
  cell_ids <- rownames(seu@meta.data[seu@meta.data[[group_var]] %in% root_cell_type, ])
  if (length(cell_ids) == 0) {
    stop("No cells found for the specified root_cell_type.")
  }
  
  cell_coords <- umap_coords[cell_ids, , drop = FALSE]
  
  # ---- 4. Extract MST vertex names and layout ----
  mst_vertex_names <- V(mst)$name
  mst_layout <- umap_coords[mst_vertex_names, , drop = FALSE]
  
  # ---- 5. Map each cell to its nearest MST vertex using FNN ----
  knn_map <- get.knnx(mst_layout, cell_coords, k = 1)
  closest_vertex_indices <- knn_map$nn.index[, 1]
  closest_vertex_names <- mst_vertex_names[closest_vertex_indices]
  
  # ---- 6. Find the most common vertex as root ----
  root_vertex <- names(which.max(table(closest_vertex_names)))
  cat("Most common root MST vertex for", root_cell_type, ":", root_vertex, "\n")
  
  
  if(isTRUE(return_root_cell)){
    return(root_vertex)
  }else{
    seu@misc$root_cell<-root_vertex
    return(seu)
  }
  
}



InferPseudotime <- function(seu, return_pseudotime = FALSE) {
  # Load necessary library
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  library(igraph)
  
  # Retrieve MST graph and root cell
  mst_graph <- seu@misc$mst_graph
  root_cell <- seu@misc$root_cell
  
  # Ensure root cell exists in the graph
  stopifnot(root_cell %in% V(mst_graph)$name)
  
  # Perform BFS traversal from root to get pseudotime distances
  bfs_result <- igraph::bfs(
    mst_graph,
    root = root_cell,
    dist = TRUE,
    order = TRUE
  )
  
  # Extract distances from BFS (pseudotime)
  pseudotime <- setNames(bfs_result$dist, V(mst_graph)$name)
  pseudotime_unnormalized <- pseudotime
  
  # Optional: remove unreachable nodes (NA distances)
  pseudotime <- pseudotime[!is.na(pseudotime)]
  
  # Normalize pseudotime to 0–1 range
  scale.factor <- max(pseudotime) - min(pseudotime)
  if (scale.factor > 0) {
    pseudotime <- pseudotime / scale.factor
  } else {
    warning("Pseudotime scale factor is 0; skipping normalization.")
  }
  
  # Return results
  if (return_pseudotime) {
    return(list(
      pseudotime = pseudotime,
      pseudotime_unnormalized = pseudotime_unnormalized
    ))
  } else {
    seu@misc$pseudotime <- pseudotime
    seu@misc$pseudotime_unnormalized <- pseudotime_unnormalized
    return(seu)
  }
}





PlotPseudotime<-function(seu,
                         reduction="umap"
){
  
  # plot the pseudotime on UMAP
  
  library(ggplot2)
  
  # Get UMAP coordinates
  umap_coords <- Embeddings(seu, reduction = reduction)
  pseudotime<-seu@misc$pseudotime
  
  
  # Combine with pseudotime
  plot_df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    pseudotime = pseudotime[rownames(umap_coords)]
  )
  
  # Remove NAs (cells without pseudotime)
  plot_df <- na.omit(plot_df)
  
  
  ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(title = "UMAP colored by Pseudotime (MST BFS)",
         color = "Pseudotime")
  
  # Add pseudotime to Seurat metadata
  seu$pseudotime <- NA
  seu$pseudotime[names(pseudotime)] <- pseudotime
  
  
  # Then plot using Seurat
  p<-FeaturePlot(seu, "pseudotime", reduction = "umap") +
    scale_color_viridis_c(option = "plasma") 
  
  p
  
  return(p)
  
}






