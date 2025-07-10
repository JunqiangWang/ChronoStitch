

gg.theme<-theme(
  # title=element_blank(),
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  #legend.position="none",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  # plot.background=element_blank(),
  plot.background = element_rect(color = "black", linewidth = 0.5)
)




gg.theme2<-theme(
  # title=element_blank(),
  axis.line=element_blank(),
  # axis.text.x=element_blank(),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  # axis.text.y=element_blank(),
  #axis.ticks=element_blank(),
  #axis.title.x=element_blank(),
  #axis.title.y=element_blank(),
  #legend.position="none",
  panel.background=element_blank(),
  #panel.border=element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
  # plot.background=element_blank(),
  # plot.background = element_rect(color = "black", linewidth = 0.1)
  
)




highlight_cells_on_umap <- function(seu, cells_by_traj, umap_reduction = "umap") {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  
  # 1. Flatten the list into a long format: (trajectory, cell)
  highlight_df <- enframe(cells_by_traj, name = "trajectory", value = "cell") %>%
    unnest(cell)
  
  # 2. Count how many trajectories each cell appears in
  cell_counts <- highlight_df %>%
    count(cell, name = "n_trajectories")
  
  # 3. Label cells that are in multiple trajectories
  highlight_df <- highlight_df %>%
    left_join(cell_counts, by = "cell") %>%
    mutate(traj_label = ifelse(n_trajectories > 1, "Multiple", trajectory))
  
  # 4. Get UMAP coordinates
  umap_coords <- Embeddings(seu, reduction = umap_reduction)
  umap_df <- data.frame(
    cell = rownames(umap_coords),
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2]
  )
  
  # 5. Merge UMAP and labels
  umap_plot_df <- left_join(umap_df, highlight_df, by = "cell")
  
  # 6. Define color palette: one for each unique trajectory label
  all_labels <- unique(umap_plot_df$traj_label)
  traj_colors <- setNames(
    c(RColorBrewer::brewer.pal(max(3, length(all_labels) - 1), "Set1"), "pink"),
    c(setdiff(all_labels, "Multiple"), "Multiple")
  )
  
  # 7. Plot
  p<-ggplot() +
    # Unassigned cells as background
    geom_point(data = umap_df[!(umap_df$cell %in% umap_plot_df$cell), ],
               aes(x = UMAP_1, y = UMAP_2), color = "grey", size = 0.1) +
    
    # Assigned cells colored by trajectory or "Multiple"
    geom_point(data = umap_plot_df,
               aes(x = UMAP_1, y = UMAP_2, color = traj_label),
               size = 0.1, alpha = 1) +
    
    theme_minimal() +
    labs(title = "",
         color = "Trajectory") +
    scale_color_manual(values = traj_colors) +
    gg.theme2+
    theme(legend.position = "none")
  
  return(p)
}




plot_trajectories_with_assignments <- function(seu, mst, trajectories, assigned_trajectory_list, umap_reduction = "umap") {
  # 1. Get UMAP coordinates
  umap_coords <- Embeddings(seu, reduction = umap_reduction)
  umap_df <- data.frame(
    cell = rownames(umap_coords),
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2]
  )
  
  # 2. Convert assigned trajectory list to long format
  assignment_df <- tibble::enframe(assigned_trajectory_list, name = "cell", value = "trajectory") %>%
    tidyr::unnest(trajectory)
  
  # 3. Identify multi-assigned cells
  multi_assignment <- assignment_df %>%
    group_by(cell) %>%
    summarize(n = n(), .groups = "drop") %>%
    mutate(multi = n > 1)
  
  assignment_df <- left_join(assignment_df, multi_assignment, by = "cell")
  assignment_df$trajectory_label <- ifelse(assignment_df$multi, "Multiple", assignment_df$trajectory)
  
  # Merge with UMAP
  plot_df <- left_join(umap_df, assignment_df, by = "cell")
  
  # 4. Extract MST paths as line segments
  path_lines_df <- do.call(rbind, lapply(seq_along(trajectories), function(i) {
    path_cells <- trajectories[[i]]
    coords <- umap_df %>%
      filter(cell %in% path_cells) %>%
      mutate(order = match(cell, path_cells),
             branch = paste0("Trajectory_", i)) %>%
      arrange(order)
    return(coords)
  }))
  
  # 5. Plot
  ggplot() +
    # Background cells not assigned
    geom_point(data = umap_df[!(umap_df$cell %in% plot_df$cell), ],
               aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.1) +
    
    # Assigned cells
    geom_point(data = plot_df,
               aes(x = UMAP_1, y = UMAP_2, color = trajectory_label), size = 0.1, alpha = 1) +
    
    # MST trajectory paths
    geom_path(data = path_lines_df,
              aes(x = UMAP_1, y = UMAP_2, group = branch),
              color = "black", size = 0.1, alpha = 1) +
    
    theme_minimal() +
    labs(title = "Cells Assigned to Trajectories with MST Overlay",
         color = "Trajectory Assignment") +
    gg.theme2
}






