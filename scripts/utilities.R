# utilities 

####################################
# MAKE PLOTTING DFS
####################################
make_df <- function(sce, dim_what, gene_name){
  
  dim1 <- as.numeric(sce@int_colData@listData[["reducedDims"]]@listData[[paste(dim_what)]][, 1])
  dim2 <- as.numeric(sce@int_colData@listData[["reducedDims"]]@listData[[paste(dim_what)]][, 2])
  dim3 <- as.numeric(sce@int_colData@listData[["reducedDims"]]@listData[[paste(dim_what)]][, 3])
  
  cluster <- sce$cluster
  celltype <- sce$label
  
  # search for gene 
  idx <- grep(paste0("\\", gene_name, "\\b"), rownames(sce), ignore.case = TRUE)
  gene <- rownames(sce)[idx]
  gene_counts <- as.vector(logcounts(sce)[idx, ])
  
  if (length(celltype) == 0) {
    df <- data.frame(dim1, dim2, dim3, cluster, gene_counts)
    
    if (dim_what == "TSNE") {
      names(df) <- c("tsne1", "tsne2", "tsne3", "cluster", paste(gene))
    } else if (dim_what == "UMAP") {
      names(df) <- c("umap1", "umap2", "umap3", "cluster", paste(gene))
    } else if (dim_what == "PCA") {
      names(df) <- c("pca1", "pca2", "pca3", "cluster", paste(gene))
    }
    
  } else {
    df <- data.frame(dim1, dim2, dim3, cluster, celltype, gene_counts)
  
    if (dim_what == "TSNE") {
      names(df) <- c("tsne1", "tsne2", "tsne3", "cluster", "celltype", paste(gene))
    } else if (dim_what == "UMAP") {
      names(df) <- c("umap1", "umap2", "umap3", "cluster", "celltype", paste(gene))
    } else if (dim_what == "PCA") {
      names(df) <- c("pca1", "pca2", "pca3", "cluster", "celltype", paste(gene))
    }}
  
  return(df)
  
}

####################################
# MAKE 2D PLOTS
####################################


make_plots_2d <- function(sce, df, color_what, dimred_type, colour_list){
  
  if (color_what == "cluster"){
    
    cols <- colour_list[1:as.numeric(length(unique(sce$cluster)))]
    
  } else if (color_what == "celltype"){
      
      cols <- colour_list[1:as.numeric(length(unique(sce$label)))]
      
  } else {cols <- "Reds"}
  
  if (dimred_type == "TSNE"){names(df)[1:3] <- c("tsne1", "tsne2", "tsne3")
  
    p <- plot_ly(data = df,
                 x = ~df$tsne1, y = ~df$tsne2,
                 opacity = 1,
                 color = ~df[, color_what],
                 type = "scatter",
                 mode = "markers",
                 marker = list(size = 3), 
                 colors = cols) %>% 
      layout(legend= list(itemsizing='constant')) } else if (dimred_type == "UMAP"){names(df)[1:3] <- c("umap1", "umap2", "umap3")
    
    p <- plot_ly(data = df,
                 x = ~df$umap1, y = ~df$umap2,
                 opacity = 1,
                 color = ~df[, color_what],
                 type = "scatter",
                 mode = "markers",
                 marker = list(size = 3), 
                 colors = cols) %>% 
      layout(legend= list(itemsizing='constant')) } else if (dimred_type == "PCA"){names(df)[1:3] <- c("pca1", "pca2", "pca3")
      
      p <- plot_ly(data = df,
                   x = ~df$pca1, y = ~df$pca2,
                   opacity = 1,
                   color = ~df[, color_what],
                   type = "scatter",
                   mode = "markers",
                   marker = list(size = 3), 
                   colors = cols) %>% 
        layout(legend= list(itemsizing='constant'))
      
      }
  
  
  return(p)
  
  }


####################################
# MAKE 3D PLOTS
####################################

make_plots_3d <- function(sce, df, color_what, dimred_type, colour_list){
  
  if (color_what == "cluster"){
    
    cols <- colour_list[1:as.numeric(length(unique(sce$cluster)))]
    
  } else if (color_what == "celltype"){
    
    cols <- colour_list[1:as.numeric(length(unique(sce$label)))]
    
  } else {cols <- "Reds"}
  
  if (dimred_type == "TSNE"){names(df)[1:3] <- c("tsne1", "tsne2", "tsne3")
  
  p <- plot_ly(data = df,
               x = ~df$tsne1, y = ~df$tsne2, z = ~df$tsne3,
               opacity = 1,
               color = ~df[, color_what],
               type = "scatter3d",
               mode = "markers",
               marker = list(size = 3), 
               colors = cols) %>% 
    layout(legend= list(itemsizing='constant')) } else if (dimred_type == "UMAP"){names(df)[1:3] <- c("umap1", "umap2", "umap3")
    
    p <- plot_ly(data = df,
                 x = ~df$umap1, y = ~df$umap2, z = ~df$umap3,
                 opacity = 1,
                 color = ~df[, color_what],
                 type = "scatter3d",
                 mode = "markers",
                 marker = list(size = 3), 
                 colors = cols) %>% 
      layout(legend= list(itemsizing='constant')) } else if (dimred_type == "PCA"){names(df)[1:3] <- c("pca1", "pca2", "pca3")
      
      p <- plot_ly(data = df,
                   x = ~df$pca1, y = ~df$pca2, z = ~df$pca3,
                   opacity = 1,
                   color = ~df[, color_what],
                   type = "scatter",
                   mode = "markers",
                   marker = list(size = 3), 
                   colors = cols) %>% 
        layout(legend= list(itemsizing='constant'))
    }
  
  
  return(p)
  
}

####################################
# MAKE BAR GRAPHS BASED ON CELL TYPE
####################################

filter_and_plot <- function(df, cutoff, table_what, gene_name){
  
  df <- df %>% 
    filter(df[, gene_name] >= cutoff)
  
  counts_df <- as.data.frame(table(df[, table_what]))
  names(counts_df) <- c("celltype", "frequency")
  
  p <- ggplot(data=counts_df, aes(x=celltype, y=frequency, fill = celltype)) +
    geom_bar(stat="identity") + 
    coord_flip() +
    scale_fill_brewer(palette="Dark2") +
    theme_bw()
  
  return(p)
  
}






































