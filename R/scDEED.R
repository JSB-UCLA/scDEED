# setup global variables 
input_counts <- nFeature_RNA <- fir_input_data <- umap_nm_table <- n <- m <- umap_n_table <- umap_m_table <- n.neighbors <- `number of dubious cells` <- tsne_table <- percent.mt <- p <- i <- NULL

# Permutation for input as Seurat object
Permuted = function (pbmc)
{
    pbmc.permuted = pbmc
    X=pbmc@assays$RNA@scale.data
    X_permuted = pbmc.permuted@assays$RNA@scale.data
    set.seed(1000)
for (i in 1:dim(X)[1])
{
    row = pracma::randperm(dim(X)[2])
    X_permuted[i,]=X[i,row]
}
pbmc.permuted@assays$RNA@scale.data = X_permuted

    return(pbmc.permuted)
}


# Find distances under PCA and tSNE for both original and permuted matrices.
Distances.PCA.big = function(pbmc,pbmc.permuted, K, perplexity_score)
{

  distances<-distances::distances
   pbmc.permuted<-Seurat::RunPCA(pbmc.permuted, npcs= K,features = Seurat::VariableFeatures(object = pbmc.permuted))


  M=pbmc@reductions[["pca"]]@cell.embeddings
  PCA_distances = distances(M[,1:K])

  M_permuted=pbmc.permuted@reductions[["pca"]]@cell.embeddings
  PCA_distances_permuted = distances(M_permuted[,1:K])


  tSNE_distances = distances(pbmc@reductions$tsne@cell.embeddings)

  pbmc.permuted <- Seurat::RunTSNE(pbmc.permuted, dims.use = 1:K, seed.use=100, perplexity = perplexity_score)
  tSNE_distances_permuted = distances(pbmc.permuted@reductions$tsne@cell.embeddings)

  results.PCA<-list("object"=pbmc,"object.permuted"=pbmc.permuted,"PCA_distances"=PCA_distances,"PCA_distances_permuted"=PCA_distances_permuted,"tSNE_distances"=tSNE_distances,"tSNE_distances_permuted"=tSNE_distances_permuted)

  return(results.PCA)
}

# Find distances under PCA and UMAP for both original and permuted matrices.
Distances.PCA.UMAP.big = function(pbmc,pbmc.permuted, K)
{

  distances<-distances::distances
  pbmc.permuted<-Seurat::RunPCA(pbmc.permuted, npcs= K,features = Seurat::VariableFeatures(object = pbmc.permuted))


  M=pbmc@reductions[["pca"]]@cell.embeddings
  PCA_distances = distances(M[,1:K])

  M_permuted=pbmc.permuted@reductions[["pca"]]@cell.embeddings
  PCA_distances_permuted = distances(M_permuted[,1:K])


  UMAP_distances = distances(pbmc@reductions$umap@cell.embeddings)

  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, dims = 1:K, seed.use=100)
  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)

  results.PCA<-list("object"=pbmc,"object.permuted"=pbmc.permuted,"PCA_distances"=PCA_distances,"PCA_distances_permuted"=PCA_distances_permuted,"UMAP_distances"=UMAP_distances,"UMAP_distances_permuted"=UMAP_distances_permuted)

  return(results.PCA)
}


# Choose optimal perplexity parameter under tSNE. Here for one parameter
ChoosePerplexity<- function(pbmc,pbmc.permuted, K, perplexity_score){
  distances<-distances::distances
   pbmc_perplex = Seurat::RunTSNE(pbmc, dims.use = 1:K, seed.use = 100, perplexity = perplexity_score, do.fast = T)
    print("tSNE done")
    tSNE_distances = distances(pbmc_perplex@reductions$tsne@cell.embeddings)
    pbmc.permuted <- Seurat::RunTSNE(pbmc.permuted, dims.use = 1:K, seed.use = 100, perplexity = perplexity_score, do.fast = T)
      print("permuted tSNE done")
    tSNE_distances_permuted = distances(pbmc.permuted@reductions$tsne@cell.embeddings)
    results.PCA<-list("object"= pbmc_perplex,"object.permuted"=pbmc.permuted,"tSNE_distances"=tSNE_distances,"tSNE_distances_permuted"=tSNE_distances_permuted)
        return(results.PCA)
}



# Choose the n.neighbors parameter in UMAP. Here for one parameter
ChoosenNeighbors<- function(pbmc,pbmc.permuted, reduction.method, K, n, m){
  distances<-distances::distances
  pbmc_nNeighbors = Seurat::RunUMAP(pbmc, reduction = reduction.method,dims = 1:K, seed.use = 100, n.neighbors=n, min.dist = m)
  print("UMAP done")
  UMAP_distances = distances(pbmc_nNeighbors@reductions$umap@cell.embeddings)
  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, reduction = reduction.method, dims = 1:K, seed.use = 100, n.neighbors=n, min.dist = m)
  print("permuted UMAP done")
  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)
  results<-list("object"= pbmc_nNeighbors,"object.permuted"=pbmc.permuted,"UMAP_distances"=UMAP_distances,"UMAP_distances_permuted"=UMAP_distances_permuted)
  return(results)
}

 

Cell.Similarity.UMAP = function(Euc_distances,Euc_distances_permuted,UMAP_distances,UMAP_distances_permuted,percent)
    {
  numberselected = floor((dim(Euc_distances)[2])*percent)
  rho_UMAP<-foreach::`%do%`(foreach::foreach(i=1:(dim(Euc_distances)[2])), cor((UMAP_distances[i,order(Euc_distances[i,])][2:(numberselected+1)]),(sort(UMAP_distances[i,])[2:(numberselected+1)])))

  rho_UMAP = as.numeric(rho_UMAP)
  print("UMAP done")
  rho_UMAP_permuted<- foreach::`%do%`(foreach::foreach(i=1:(dim(Euc_distances)[2]) ), cor((UMAP_distances_permuted[i,order(Euc_distances_permuted[i,])][2:(numberselected+1)]),(sort(UMAP_distances_permuted[i,])[2:(numberselected+1)])))
  rho_UMAP_permuted = as.numeric(rho_UMAP_permuted)
  print("permuted UMAP done")
  similarity_results_UMAP <-list("rho_UMAP"=rho_UMAP,"rho_UMAP_permuted"=rho_UMAP_permuted)
  return(similarity_results_UMAP)
}


Cell.Similarity.tSNE = function(Euc_distances,Euc_distances_permuted,tSNE_distances,tSNE_distances_permuted,percent)
{
    numberselected = floor((dim(Euc_distances)[2])*percent)

    rho_tSNE<- foreach::`%do%`(foreach::foreach(i=1:(dim(Euc_distances)[2])), cor((tSNE_distances[i,order(Euc_distances[i,])][2:(numberselected+1)]),(sort(tSNE_distances[i,])[2:(numberselected+1)])))
    rho_tSNE = as.numeric(rho_tSNE)
    print("tSNE done")
    rho_tSNE_permuted<- foreach::`%do%`(foreach::foreach(i=1:(dim(Euc_distances)[2])), cor((tSNE_distances_permuted[i,order(Euc_distances_permuted[i,])][2:(numberselected+1)]),(sort(tSNE_distances_permuted[i,])[2:(numberselected+1)])))
    rho_tSNE_permuted = as.numeric(rho_tSNE_permuted)
    print("permuted tSNE done")
       similarity_results_tSNE <-list("rho_tSNE"=rho_tSNE,"rho_tSNE_permuted"=rho_tSNE_permuted)
    return(similarity_results_tSNE)
}



Cell.Classify.tSNE = function(rho_tSNE,rho_tSNE_permuted){
    y = seq(1,length(rho_tSNE),length=length(rho_tSNE))

    rho_tSNE_upper = quantile(rho_tSNE_permuted,0.95)
    rho_tSNE_lower = quantile(rho_tSNE_permuted,0.05)


    tSNE_badindex = which(rho_tSNE < rho_tSNE_lower)
    tSNE_intindex = which(rho_tSNE< rho_tSNE_upper&rho_tSNE >= rho_tSNE_lower)
    tSNE_goodindex = setdiff(y,union(tSNE_badindex,tSNE_intindex))
    ClassifiedCells_tSNE.results<-list("tSNE_badindex"=tSNE_badindex,"tSNE_intindex"=tSNE_intindex,"tSNE_goodindex"=tSNE_goodindex)
    return(ClassifiedCells_tSNE.results)
}

Cell.Classify.UMAP = function(rho_UMAP,rho_UMAP_permuted){
    y = seq(1,length(rho_UMAP),length=length(rho_UMAP))
    rho_UMAP_upper = quantile(rho_UMAP_permuted,0.95)
    rho_UMAP_lower = quantile(rho_UMAP_permuted,0.05)

    UMAP_badindex = which(rho_UMAP < rho_UMAP_lower)
    UMAP_intindex = which(rho_UMAP< rho_UMAP_upper&rho_UMAP >= rho_UMAP_lower)
    UMAP_goodindex = setdiff(y,union(UMAP_badindex,UMAP_intindex))

   ClassifiedCells_UMAP.results<-list("UMAP_badindex"=UMAP_badindex,"UMAP_intindex"=UMAP_intindex,"UMAP_goodindex"=UMAP_goodindex)
    return(ClassifiedCells_UMAP.results)
}

chooseK = function(input_data){
  if(class(input_data) != "Seurat"){
    dupli <- duplicated(colnames(input_data))
    if(length(which(dupli == TRUE)) != 0){
      input_data <- input_data[, !dupli]
    }
    negative <- which(as.matrix(input_data < 0), arr.ind = T)[, 2]
    if(length(negative) != 0){
      warning("There are negative values in the count matrix, and we have removed the corresponding cells and processd.")
      input_data <- input_data[, -negative]
    }
    
    zero_row <- which(Matrix::rowSums(input_data) == 0)
    zero_col <- which(Matrix::colSums(input_data) == 0)
    
    
    
    if(length(which(zero_row == TRUE)) != 0){
      warning("There are rows containing all zeros, and we have removed them correspondingly and proceed.")
      input_data <- input_data[-zero_row, ]
    }
    if(length(which(zero_col == TRUE)) != 0){
      warning("There are columns containing all zeros, and we have removed them correspondingly and proceed.")
      input_data <- input_data[, -zero_col]
    }
    input_data = suppressWarnings({Seurat::CreateSeuratObject(counts = input_data)})
  }
  input_data <- Seurat::FindVariableFeatures(input_data)
  input_data <- Seurat::ScaleData(input_data)
  input_data <- Seurat::RunPCA(input_data, features = Seurat::VariableFeatures(object = input_data))
  print(Seurat::ElbowPlot(input_data))
}

# helper function to run ChoosenNeighbors, Cell.Similarity.UMAP, and Cell.Classify.UMAP in umap_tsne_process
umap_optimize = function(input_data,input_data.permuted, reduction.method, K, n, m, results.PCA, similarity_percent){
  results<-ChoosenNeighbors(pbmc = input_data ,pbmc.permuted = input_data.permuted, reduction.method = reduction.method, K = K, n = n, m = m)
  
  similarity_score_UMAP <-Cell.Similarity.UMAP(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,results$UMAP_distances,results$UMAP_distances_permuted, similarity_percent)

  ClassifiedCells_UMAP<-Cell.Classify.UMAP(similarity_score_UMAP$rho_UMAP,similarity_score_UMAP$rho_UMAP_permuted)
  
  return(length(ClassifiedCells_UMAP$UMAP_badindex))
}

tsne_optimize = function(input_data, input_data.permuted, num_pc, perplexity, results.PCA, similarity_percent){
  results<-ChoosePerplexity(input_data,input_data.permuted, num_pc, perplexity)
    similarity_score_tSNE <-Cell.Similarity.tSNE(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,results$tSNE_distances,results$tSNE_distances_permuted, similarity_percent)
    ClassifiedCells_tSNE<-Cell.Classify.tSNE(similarity_score_tSNE$rho_tSNE,similarity_score_tSNE$rho_tSNE_permuted)
   return(length(ClassifiedCells_tSNE$tSNE_badindex))
}

scDEED = function(input_data, num_pc, n_neighbors = c(seq(from=5,to=30,by=1),35,40,45,50), min.dist = seq(0.1,0.9, by = 0.2), similarity_percent = 0.5, visualization = FALSE, use_method = "umap",
                             perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)), perplexity_score = 30, optimize_neib = TRUE ,optimize_min = TRUE){
  # check if user just want to see the demo output 
  

  demo_tsne <- (identical(input_data,input_counts) && perplexity == c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)))
  demo_umap_mn <- (identical(input_data,input_counts) && n_neighbors == c(seq(from=5,to=30,by=1),35,40,45,50) && min.dist == seq(0.1,0.9, by = 0.2))
  demo_umap_n <- (identical(input_data,input_counts) && n_neighbors == c(seq(from=5,to=30,by=1),35,40,45,50))
  demo_umap_m <- (identical(input_data,input_counts) && min.dist == seq(0.1,0.9, by = 0.2))
  
  if(class(input_data) != "Seurat"){
    dupli <- duplicated(colnames(input_data))
    if(length(which(dupli == TRUE)) != 0){
      input_data <- input_data[, !dupli]
    }
    negative <- which(as.matrix(input_data < 0), arr.ind = T)[, 2]
    if(length(negative) != 0){
      warning("There are negative values in the count matrix, and we have removed the corresponding cells and processd.")
      input_data <- input_data[, -negative]
    }
    
    zero_row <- which(Matrix::rowSums(input_data) == 0)
    zero_col <- which(Matrix::colSums(input_data) == 0)
    
     
    
    if(length(which(zero_row == TRUE)) != 0){
      warning("There are rows containing all zeros, and we have removed them correspondingly and proceed.")
      input_data <- input_data[-zero_row, ]
    }
    if(length(which(zero_col == TRUE)) != 0){
      warning("There are columns containing all zeros, and we have removed them correspondingly and proceed.")
      input_data <- input_data[, -zero_col]
    }
    
  suppressWarnings({input_data = Seurat::CreateSeuratObject(counts = input_data)})
  }
  
  
  
  input_data[["percent.mt"]] <- Seurat::PercentageFeatureSet(input_data, pattern = "^MT-")
  input_data <- subset(input_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  input_data <- Seurat::NormalizeData(input_data, normalization.method = "LogNormalize", scale.factor = 10000)
  input_data <- Seurat::NormalizeData(input_data)
  
  input_data <- Seurat::FindVariableFeatures(input_data)

  input_data <- Seurat::ScaleData(input_data)

  input_data <- Seurat::RunPCA(input_data, features = Seurat::VariableFeatures(object = input_data))
  nsamples <- ncol(input_data)
  
  # subsampling process
  if(nrow(input_data[[]]) > 30000){
    full_cell <- Seurat::Cells(input_data)
    if(length(levels(input_data)) <= 1){
      # no cell types 
      sub_cells <- sample(full_cell, 30000, replace = FALSE)
      input_data <- subset(input_data, cells = sub_cells)
    }
    else{
      # with cell types
      cellTypes <- levels(input_data)
      preserved_cells_num <- (100/length(cellTypes))*0.01*length(full_cell)
      # sampling for each cell type
      
      # start with the first cell type
      input_data <- subset(input_data, idents = cellTypes[1])
      fir_cells <- Seurat::Cells(fir_input_data)
      if(length(fir_cells) > preserved_cells_num){
        sub_cells <- sample(fir_cells, preserved_cells_num, replace = FALSE)
        input_data <- subset(input_data, cells = sub_cells)
      }
      for (i in cellTypes[-1]) {
        next_input_data <- subset(input_data, idents = i)
        next_cells <- Seurat::Cells(next_input_data)
        if(length(next_cells) > preserved_cells_num){
          sub_cells <- sample(next_cells, preserved_cells_num, replace = FALSE)
          next_input_data <- subset(next_input_data, cells = sub_cells)
        }
        input_data <- merge(input_data, next_input_data)
      }
    }
    # add the features and pca to input_data
    input_data <- Seurat::FindVariableFeatures(input_data)
    
    input_data <- Seurat::ScaleData(input_data)
    
    input_data <- Seurat::RunPCA(input_data, features = Seurat::VariableFeatures(object = input_data))
  }
  
  if(use_method == "umap"){
    
   
    if (is.null(input_data@reductions$umap)){
      input_data <- Seurat::RunUMAP(input_data, dims = 1:num_pc)
    }
    input_data.permuted<-Permuted(input_data)
    input_data.permuted<-Seurat::RunPCA(input_data.permuted, npcs= num_pc, features = Seurat::VariableFeatures(object = input_data.permuted))

    results.PCA <- Distances.PCA.UMAP.big(input_data, input_data.permuted, K=num_pc)

    
    # set up the default parameters for runumap
    default_neib <- 30L
    default_min <- 0.3
    
    # set up the parameters used in the final runumap.
    final_neib <- default_neib
    final_min <- default_min
    
    if(optimize_neib == TRUE && optimize_min == TRUE){
      # quick output umap sample result if user input the sample data
      if(demo_umap_mn){

        dubious_number_UMAP <- umap_nm_table
        dub_para <- umap_nm_table
        all_dub <- umap_nm_table$`number of dubious cells`
      }
      else{
        all_pairs <- expand.grid(n_neighbors, min.dist)
        
        all_dub <-foreach::`%do%`(foreach::foreach(n = all_pairs$Var1, m = all_pairs$Var2, .combine = "c"), umap_optimize(input_data = input_data, input_data.permuted = input_data.permuted,
                                                                                                                          reduction.method = "pca", K = num_pc,
                                                                                                                          n = n, m = m, results.PCA = results.PCA, similarity_percent = similarity_percent) 
        )
        dubious_number_UMAP <- cbind(all_pairs, all_dub)
        
        dub_para <- data.frame("n.neighbors" = all_pairs$Var1, "min.dist" = all_pairs$Var2, "number of dubious cells" = all_dub)
        colnames(dub_para) <- c("n.neighbors", "min.dist", "number of dubious cells")
      }
 
    
      
      
      best_para <- dubious_number_UMAP[which(all_dub == min(all_dub)),c(1, 2)]
      if(!is.null(nrow(best_para) )){
        row_sum <- rowSums(best_para)
        best_para <- best_para[which(row_sum == min(row_sum)) , ]
      }
      colnames(best_para) <- c("n.neighbers", "min.dist")

      final_neib <- best_para$n.neighbers
      final_min <- best_para$min.dist
    }
    
    if(optimize_neib == TRUE && optimize_min == FALSE){
      
      # quick result
      if(demo_umap_n){
    
        dubious_number_UMAP_neib <- umap_n_table$`number of dubious cells`
      }
      else{
        dubious_number_UMAP_neib <-foreach::`%do%`(foreach::foreach(n = n_neighbors, .combine = "c"), umap_optimize(input_data = input_data, input_data.permuted = input_data.permuted,
                                                                                                                    reduction.method = "pca", K = num_pc,
                                                                                                                    n = n, m = default_min, results.PCA = results.PCA, similarity_percent = similarity_percent) 
        )
      }
      
      
      best_para_neib <- n_neighbors[which(dubious_number_UMAP_neib == min(dubious_number_UMAP_neib))]
      
      
      if(length(best_para_neib) != 0){
        best_para_neib <- min(best_para_neib)
      }
      dub_neighbor <- data.frame("n.neighbors" = n_neighbors, "number of dubious cells" = dubious_number_UMAP_neib)
      colnames(dub_neighbor) <- c("n.neighbors", "number of dubious cells")
      final_neib <- best_para_neib
    }
    
    
    if(optimize_min == TRUE && optimize_neib == FALSE){
    
      # quick result
      if(demo_umap_m){
   
        dubious_number_UMAP_min <- umap_m_table$`number of dubious cells`
      }
      else{
        dubious_number_UMAP_min <- foreach::`%do%`(foreach::foreach(m = min.dist, .combine = "c"), umap_optimize(input_data = input_data, input_data.permuted = input_data.permuted,
                                                                                                                 reduction.method = "pca", K = num_pc,
                                                                                                                 n = default_neib, m = m, results.PCA = results.PCA, similarity_percent = similarity_percent) 
        )
      }
      
      
      best_para_min <- min.dist[which(dubious_number_UMAP_min == min(dubious_number_UMAP_min))]
      if(length(best_para_min) != 0){
        best_para_min <- min(best_para_min)
      }
      
      dub_min_dist <- data.frame("min.dist" = min.dist, "number of dubious cells" = dubious_number_UMAP_min)
      colnames(dub_min_dist) <- c("min.dist", "number of dubious cells")
      final_min <- best_para_min
    }
    
    
  
    res <- ChoosenNeighbors(input_data, input_data.permuted, "pca", num_pc, n = final_neib, m = final_min)
    similarity_score_UMAP <-Cell.Similarity.UMAP(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,res$UMAP_distances,res$UMAP_distances_permuted, similarity_percent)
    
    ClassifiedCells_UMAP<-Cell.Classify.UMAP(similarity_score_UMAP$rho_UMAP,similarity_score_UMAP$rho_UMAP_permuted)
    cell_list <- rownames(res$object@meta.data)
    if(visualization == TRUE){
      
      bad_graph <- Seurat::DimPlot(res$object, reduction = "umap", cells.highlight = list(`Dubious cells`= ClassifiedCells_UMAP$UMAP_badindex), cols.highlight = "red") 
      levels(bad_graph$data$highlight)[match("Unselected",levels(bad_graph$data$highlight))] <- "Other Cells"
      trust_graph <- Seurat::DimPlot(res$object, reduction = "umap", cells.highlight = list(`Trustworthy cells`= ClassifiedCells_UMAP$UMAP_goodindex), cols.highlight = "blue") 
      levels(trust_graph$data$highlight)[match("Unselected",levels(trust_graph$data$highlight))] <- "Other Cells"
      

      if(optimize_neib == TRUE && optimize_min == FALSE){
        highlight_neib <- subset(dub_neighbor, n.neighbors == best_para_neib)
        input_data_dubious_plot_neib <- ggplot2::ggplot(data=dub_neighbor, ggplot2::aes(x=n.neighbors, y=`number of dubious cells`, group=1))+
          ggplot2::geom_point(size = 5)+
          ggplot2::geom_point(data=highlight_neib, ggplot2::aes(x=n.neighbors, y=`number of dubious cells`), color='cyan',size = 5) +
          ggplot2::geom_vline(xintercept=highlight_neib$n.neighbors, linetype="dotted") +
          ggplot2::annotate(geom = "text", x =highlight_neib$n.neighbors, y = max(dub_neighbor$`number of dubious cells`), label = "optimized",color='cyan',size =4, vjust = "inward", hjust = "inward")+
          ggplot2::labs(x = "n.neighbors", y = "# of dubious cell embeddings") + ggplot2::theme_bw() +
          ggplot2::theme(text=ggplot2::element_text(size=20),panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),axis.text=ggplot2::element_text(size=20), axis.line = ggplot2::element_line(colour = "black"))
        output <-list(dub_neighbor, best_para_neib, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex], bad_graph, trust_graph, input_data_dubious_plot_neib)
        
        names(output) <- c("number of dubious cells corresponding to n.neighbors list", "best n.neighbors", 
                           "list of dubious cells corresponding to best n.neighbors",
                           "list of trustworthy cells corresponding to best n.neighbors",
                           "UMAP plot with dubious cells - best n.neighbors", 
                           "UMAP plot with trustworthy cells - best n.neighbors", "plot. # of dubious embeddings vs parameters")
      }
       else if(optimize_neib == FALSE && optimize_min == TRUE){
        highlight_min <- subset(dub_min_dist, min.dist == best_para_min)
        input_data_dubious_plot_min <- ggplot2::ggplot(data=dub_min_dist, ggplot2::aes(x=min.dist, y=`number of dubious cells`, group=1))+
          ggplot2::geom_point(size = 5)+
          ggplot2::geom_point(data=highlight_min, ggplot2::aes(x=min.dist, y=`number of dubious cells`), color='cyan',size = 5) +
          ggplot2::geom_vline(xintercept=highlight_min$min.dist, linetype="dotted") +
          ggplot2::annotate(geom = "text", x =highlight_min$min.dist, y = max(dub_min_dist$`number of dubious cells`), label = "optimized",color='cyan',size =4, vjust = "inward", hjust = "inward")+
          ggplot2::labs(x = "min.dist", y = "# of dubious cell embeddings") + ggplot2::theme_bw() +
          ggplot2::theme(text=ggplot2::element_text(size=20),panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),axis.text=ggplot2::element_text(size=20), axis.line = ggplot2::element_line(colour = "black"))
        output <-list(dub_min_dist, best_para_min, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex], bad_graph, trust_graph, input_data_dubious_plot_min)
        
        names(output) <- c("number of dubious cells corresponding to min.dist list", "best min.dist", 
                           "list of dubious cells corresponding to best min.dist",
                           "list of trustworthy cells corresponding to best min.dist",
                           "UMAP plot with dubious cells - best min.dist", 
                           "UMAP plot with trustworthy cells - best min.dist", "plot. # of dubious embeddings vs parameters")
      }
      else if(optimize_neib == TRUE && optimize_min == TRUE){
        highlight <- subset(dub_para, n.neighbors == best_para[1, 1] & min.dist == best_para[1, 2])
        input_data_dubious_plot <- ggplot2::ggplot(data = dub_para, ggplot2::aes(x = factor(n.neighbors), y = factor(min.dist), fill = `number of dubious cells`))  + ggplot2::geom_tile() +

          ggplot2::scale_fill_gradient(low="blue", high="red") +
          ggplot2::geom_vline(xintercept = factor(highlight$n.neighbors)) + ggplot2::annotate(geom = "text", x = factor(highlight$n.neighbors), y = factor(max(dub_para[, 2])), label = "optimized",color='cyan',size =4, vjust = "inward", hjust = "inward") +
          ggplot2::geom_hline(yintercept = factor(highlight$min.dist)) + ggplot2::annotate(geom = "text", x = factor(max(dub_para[, 1])), y = factor(highlight$min.dist), label = "optimized",color='cyan',size =4, vjust = "inward", hjust = "inward") +
          ggplot2::labs(x = "n.neighbors", y = "min.dist") + ggplot2::theme_bw() +
          ggplot2::theme(text=ggplot2::element_text(size=20),panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),axis.text=ggplot2::element_text(size=5), axis.line = ggplot2::element_line(colour = "black"))
        
        output <-list(dub_para, best_para, 
                      cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex], 
                      bad_graph, trust_graph, 
          input_data_dubious_plot)
        
        names(output) <- c("number of dubious cells corresponding to pair of n.neighbors and min.dist list", "best pair of n.neighbors and min.dist",
                           "list of dubious cells corresponding to best pair of n.neighbors and min.dist",
                           "list of trustworthy cells corresponding to best pair of n.neighbors and min.dist",
                           "UMAP plot with dubious cells - best pair of n.neighbors and min.dist", 
                           "UMAP plot with trustworthy cells - best pair of n.neighbors and min.dist", 
          "plot. # of dubious embeddings vs pair of n.neighbors and min.dist")
      }
      else{
        output <-list(cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex], 
                      bad_graph, trust_graph)
        names(output) <- c("list of dubious cells",
                           "list of trustworthy cells",
                           "UMAP plot with dubious cells", 
                           "UMAP plot with trustworthy cells")
      }
      
      
      return(output)
       }
    else{
      if(optimize_neib == TRUE && optimize_min == FALSE){
      output <-list(dub_neighbor, best_para_neib, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
      names(output) <- c("number of dubious cells corresponding to n.neighbors list", "best n.neighbors", 
                         "list of dubious cells corresponding to best n.neighbors",
                         "list of trustworthy cells corresponding to best n.neighbors")}
      else if(optimize_neib == FALSE && optimize_min == TRUE){
        output <-list(dub_min_dist, best_para_min, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
        names(output) <- c("number of dubious cells corresponding to min.dist list", "best min.dist", 
                           "list of dubious cells corresponding to best min.dist",
                           "list of trustworthy cells corresponding to best min.dist")}
      else if(optimize_neib == TRUE && optimize_min == TRUE){
        output <-list(dub_para, best_para, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
        names(output) <- c("number of dubious cells corresponding to pair of n.neighbors and min.dist list", 
                           "best pair of n.neighbors and min.dist",
                           "list of dubious cells corresponding to best pair of n.neighber and min.dist",
                           "list of trustworthy cells corresponding to best pair of n.neighber and min.dist")}
      else{
        output <-list(cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
        names(output) <- c(
                           "list of dubious cells corresponding to best n.neighber and min.dist",
                           "list of trustworthy cells corresponding to best n.neighber and min.dist")}
      return(output)
      }
      
      }

    

  
  if(use_method == "tsne"){
    if(nsamples < 91){
      stop("sample size too small")
    }
    
    if(is.null(input_data@reductions$tsne)){
      input_data <- Seurat::RunTSNE(input_data)
    }

    input_data.permuted <- Permuted(input_data)
    results.PCA <- Distances.PCA.big(input_data, input_data.permuted, K = num_pc, perplexity_score = perplexity_score)

    input_data.permuted<-Seurat::RunPCA(input_data.permuted, npcs= num_pc,features = Seurat::VariableFeatures(object = input_data.permuted))
    perplexity <- sort(perplexity)
    perplexity <- perplexity[perplexity <= floor((nsamples -1)/3)]
    

    

    
    if(demo_tsne){

      dubious_number_tSNE <- tsne_table$`number of dubious cells`
    }
    else{
      dubious_number_tSNE <-foreach::`%do%`(foreach::foreach(p = perplexity, .combine = "c"), tsne_optimize(input_data, input_data.permuted, num_pc, perplexity = p, results.PCA, similarity_percent) 
      )
    }
    
    best_para <- perplexity[which(dubious_number_tSNE == min(dubious_number_tSNE))]
    if(length(best_para) != 0){
      best_para <- min(best_para)
    }
    dub_perplex <- data.frame("perplexity" = perplexity, "number of dubious cells" = dubious_number_tSNE)
    colnames(dub_perplex) <- c("perplexity", "number of dubious cells")
    res <- ChoosePerplexity(input_data,input_data.permuted, num_pc, best_para)
    similarity_score_tSNE <-Cell.Similarity.tSNE(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,res$tSNE_distances,res$tSNE_distances_permuted,similarity_percent)
    
    ClassifiedCells_tSNE<-Cell.Classify.tSNE(similarity_score_tSNE$rho_tSNE,similarity_score_tSNE$rho_tSNE_permuted)
    cell_list <- rownames(res$object@meta.data)
    if(visualization == TRUE){
      
      bad_graph <- Seurat::DimPlot(res$object, reduction = "tsne", cells.highlight = list(`Dubious cells`= ClassifiedCells_tSNE$tSNE_badindex), cols.highlight = "red") 
      levels(bad_graph$data$highlight)[match("Unselected",levels(bad_graph$data$highlight))] <- "Other Cells"
      trust_graph <- Seurat::DimPlot(res$object, reduction = "tsne", cells.highlight = list(`Trustworthy cells`= ClassifiedCells_tSNE$tSNE_goodindex), cols.highlight = "blue") 
      levels(trust_graph$data$highlight)[match("Unselected",levels(trust_graph$data$highlight))] <- "Other Cells"
      
      # input_data_dubious <- data.frame(perplexity, dubious_number_tSNE)
      highlight <- subset(dub_perplex, perplexity == best_para)
      input_data_dubious_plot <- ggplot2::ggplot(data=dub_perplex, ggplot2::aes(x=perplexity, y=`number of dubious cells`, group=1))+
        ggplot2::geom_point(size = 5)+
        ggplot2::geom_point(data=highlight, ggplot2::aes(x=perplexity, y=`number of dubious cells`), color='cyan',size = 5) +
        ggplot2::geom_vline(xintercept=highlight$perplexity, linetype="dotted") +
        ggplot2::annotate(geom = "text", x =highlight$perplexity, y = max(dub_perplex$`number of dubious cells`), label = "optimized",color='cyan',size =4, vjust = "inward", hjust = "inward")+
        ggplot2::labs(x = "perplexity", y = "# of dubious cell embeddings") + ggplot2::theme_bw() +
        ggplot2::theme(text=ggplot2::element_text(size=20),panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),axis.text=ggplot2::element_text(size=20), axis.line = ggplot2::element_line(colour = "black"))
      
      output <-list(dub_perplex, best_para, cell_list[ClassifiedCells_tSNE$tSNE_badindex], cell_list[ClassifiedCells_tSNE$tSNE_goodindex],bad_graph, trust_graph, input_data_dubious_plot)
      names(output) <- c("number of dubious cells corresponding to perplexity list", "best perplexity", 
                         "list of dubious cell corresponding to best perplexity",
                         "list of trustworthy cell corresponding to best perplexity",
                         "tSNE plot with dubious cells - best perplexity", 
                         "tSNE plot with trustworthy cells - best perplexity",
                         "plot. # of dubious embeddings vs parameters")
      return(output)
    }
      else{
        output <-list(dub_perplex, best_para, cell_list[ClassifiedCells_tSNE$tSNE_badindex], cell_list[ClassifiedCells_tSNE$tSNE_goodindex])
        names(output) <- c("number of dubious cells corresponding to perplexity list", "best perplexity",
                           "list of dubious cell corresponding to best perplexity",
                           "list of trustworthy cell corresponding to best perplexity")
        return(output)
      }

  }

}





