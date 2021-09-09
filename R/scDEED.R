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


Distances.PCA = function(pbmc,pbmc.permuted, K)
{

pbmc <-Seurat::RunPCA(pbmc, npcs= K, features = Seurat::VariableFeatures(object = pbmc))
pbmc.permuted<-Seurat::RunPCA(pbmc.permuted, npcs= K,features = Seurat::VariableFeatures(object = pbmc.permuted))


M=pbmc@reductions[["pca"]]@cell.embeddings
PCA_distances = as.matrix(dist(M[,1:K]))

M_permuted=pbmc.permuted@reductions[["pca"]]@cell.embeddings
PCA_distances_permuted = as.matrix(dist(M_permuted[,1:K]))

results.PCA<-list("object"=pbmc,"object.permuted"=pbmc.permuted,"PCA_distances"=PCA_distances,"PCA_distances_permuted"=PCA_distances_permuted)

return(results.PCA)
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
ChoosenNeighbors<- function(pbmc,pbmc.permuted, reduction.method, K, n){
  distances<-distances::distances
  pbmc_nNeighbors = Seurat::RunUMAP(pbmc, reduction = reduction.method,dims = 1:K, seed.use = 100, n.neighbors=n, do.fast = T)
  print("UMAP done")
  UMAP_distances = distances(pbmc_nNeighbors@reductions$umap@cell.embeddings)
  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, reduction = reduction.method, dims = 1:K, seed.use = 100, n.neighbors=n, do.fast = T)
  print("permuted UMAP done")
  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)
  results<-list("object"= pbmc_nNeighbors,"object.permuted"=pbmc.permuted,"UMAP_distances"=UMAP_distances,"UMAP_distances_permuted"=UMAP_distances_permuted)
  return(results)
}

globalVariables("i")

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

chooseK = function(pbmc){
  if(class(pbmc) != "Seurat"){
    pbmc = suppressWarnings({Seurat::CreateSeuratObject(counts = pbmc)})
  }
  pbmc <- Seurat::FindVariableFeatures(pbmc)
  pbmc <- Seurat::ScaleData(pbmc)
  pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
  print(Seurat::ElbowPlot(pbmc))
}

umap_tsne_process = function(pbmc, num_pc, n_neighbors = c(seq(from=5,to=30,by=1),35,40,45,50),  similarity_percent = 0.5, visualization = FALSE, use_method = "umap",
                             perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)),perplexity_score = 30){
  if(class(pbmc) != "Seurat"){
  suppressWarnings({pbmc = Seurat::CreateSeuratObject(counts = pbmc)})
  }
  pbmc <- Seurat::FindVariableFeatures(pbmc)

  pbmc <- Seurat::ScaleData(pbmc)

  pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))

  if(use_method == "umap"){
    if (is.null(pbmc@reductions$umap)){
      pbmc <- Seurat::RunUMAP(pbmc, dims = 1:num_pc)
    }
    pbmc.permuted<-Permuted(pbmc)
    pbmc.permuted<-Seurat::RunPCA(pbmc.permuted, npcs= num_pc, features = Seurat::VariableFeatures(object = pbmc.permuted))

    results.PCA <- Distances.PCA.UMAP.big(pbmc, pbmc.permuted, K=num_pc)

    dubious_number_UMAP = rep(0,length(n_neighbors))

    for (i in 1:length(n_neighbors)){
       results<-ChoosenNeighbors(pbmc,pbmc.permuted, "pca", num_pc, n_neighbors[i])

      similarity_score_UMAP <-Cell.Similarity.UMAP(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,results$UMAP_distances,results$UMAP_distances_permuted, similarity_percent)

      ClassifiedCells_UMAP<-Cell.Classify.UMAP(similarity_score_UMAP$rho_UMAP,similarity_score_UMAP$rho_UMAP_permuted)

      dubious_number_UMAP[i] = length(ClassifiedCells_UMAP$UMAP_badindex)

    }
    best_para <- n_neighbors[which(dubious_number_UMAP == min(dubious_number_UMAP))]

    if(visualization == TRUE){
      res <- ChoosenNeighbors(pbmc, pbmc.permuted, "pca", num_pc, best_para)
      similarity_score_UMAP <-Cell.Similarity.UMAP(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,res$UMAP_distances,res$UMAP_distances_permuted, similarity_percent)

      ClassifiedCells_UMAP<-Cell.Classify.UMAP(similarity_score_UMAP$rho_UMAP,similarity_score_UMAP$rho_UMAP_permuted)
      bad_graph <- Seurat::DimPlot(res$object, reduction = "umap", cells.highlight = list(Dubious = ClassifiedCells_UMAP$UMAP_badindex)) + ggplot2::scale_color_manual(labels = c("Other cells", "Dubious cells"), values = c("grey", "red"))

      trust_graph <- Seurat::DimPlot(res$object, reduction = "umap", cells.highlight = list(Trustworthy = ClassifiedCells_UMAP$UMAP_goodindex)) + ggplot2::scale_color_manual(labels = c("Other cells", "Trustworthy cells"), values = c("grey", "blue"))
      output <-list(dubious_number_UMAP, best_para, bad_graph, trust_graph)
      names(output) <- c("dubious numbers corresponding to n.neighbors list", "best n.neighbors", "UMAP plot with dubious cells", "UMAP plot with trustworthy cells")
      return(output)
       }
    else{
      output <-list(dubious_number_UMAP, best_para)
      names(output) <- c("dubious numbers corresponding to n.neighbors list", "best n.neighbors")
      return(output)
    }

  }
  if(use_method == "tsne"){
    if(is.null(pbmc@reductions$tsne)){
      pbmc <- Seurat::RunTSNE(pbmc)
    }

    pbmc.permuted <- Permuted(pbmc)
    results.PCA <- Distances.PCA.big(pbmc, pbmc.permuted, K = num_pc, perplexity_score = perplexity_score)
    dubious_number_tSNE = rep(0,length(perplexity))
    pbmc.permuted<-Seurat::RunPCA(pbmc.permuted, npcs= num_pc,features = Seurat::VariableFeatures(object = pbmc.permuted))
    for (i in 1:length(perplexity)){
      results<-ChoosePerplexity(pbmc,pbmc.permuted, num_pc, perplexity[i])
      similarity_score_tSNE <-Cell.Similarity.tSNE(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,results$tSNE_distances,results$tSNE_distances_permuted, similarity_percent)
      ClassifiedCells_tSNE<-Cell.Classify.tSNE(similarity_score_tSNE$rho_tSNE,similarity_score_tSNE$rho_tSNE_permuted)
      dubious_number_tSNE[i] = length(ClassifiedCells_tSNE$tSNE_badindex)
    }
    best_para <- perplexity[which(dubious_number_tSNE == min(dubious_number_tSNE))]
    if(visualization == TRUE){
      res <- ChoosePerplexity(pbmc,pbmc.permuted, num_pc, best_para)
      similarity_score_tSNE <-Cell.Similarity.tSNE(results.PCA$PCA_distances,results.PCA$PCA_distances_permuted,res$tSNE_distances,res$tSNE_distances_permuted,similarity_percent)

      ClassifiedCells_tSNE<-Cell.Classify.tSNE(similarity_score_tSNE$rho_tSNE,similarity_score_tSNE$rho_tSNE_permuted)
      bad_graph <- Seurat::DimPlot(res$object, reduction = "tsne", cells.highlight = list(dubious = ClassifiedCells_tSNE$tSNE_badindex)) + ggplot2::scale_color_manual(labels = c("Other cells", "Dubious cells"), values = c("grey", "red"))
      trust_graph <- Seurat::DimPlot(res$object, reduction = "tsne", cells.highlight = list(trustworthy = ClassifiedCells_tSNE$tSNE_goodindex)) + ggplot2::scale_color_manual(labels = c("Other cells", "Trustworthy cells"), values = c("grey", "blue"))
      output <-list(dubious_number_tSNE, best_para, bad_graph, trust_graph)
      names(output) <- c("dubious numbers corresponding to perplexities", "best perplexity", "tSNE plot with dubious cells", "tSNE plot with trustworthy cells")
      return(output)
    }
      else{
        output <-list(dubious_number_tSNE, best_para)
        names(output) <- c("dubious numbers corresponding to perplexities", "best perplexity")
        return(output)
      }

  }

}



