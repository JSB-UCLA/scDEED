
# Permutation for input as Seurat object
Permuted = function (pbmc, default_assay = 'active.assay', slot = 'scale.data', K)
{
  if(default_assay == 'active.assay'){
    default_assay = pbmc@active.assay
    }
  pb <- txtProgressBar(min = 0, max = 1, initial = 0, char = "=", 
                       width = NA, title, label, style = 3, file = "")
  setTxtProgressBar(pb, 0.1)
  
  Seurat::DefaultAssay(object = pbmc) <- default_assay
  pbmc.permuted = pbmc
  X             = Seurat::GetAssayData(pbmc, slot = slot)
  X_permuted    = Seurat::GetAssayData(pbmc.permuted, slot = slot)
  
  setTxtProgressBar(pb, 0.4) 
  
  set.seed(1000)
  pb_add  <- 1 / (dim(X)[1]/0.4)
  curr_pb <- 0.4
  
  for (i in 1:dim(X)[1])
  {
    row = pracma::randperm(dim(X)[2])
    X_permuted[i,]=X[i,row]
    
    curr_pb <- curr_pb + pb_add
    setTxtProgressBar(pb, curr_pb)
  }
  
  pbmc.permuted = Seurat::SetAssayData(pbmc, slot = slot, new.data = X_permuted, assay = default_assay)
  pbmc.permuted <- Seurat::RunPCA(pbmc.permuted, npcs= K,features = Seurat::VariableFeatures(object = pbmc.permuted))
  
  setTxtProgressBar(pb, 1)
  return(pbmc.permuted)
}

#this function calculates the distances in the pre-embedding space for the original and permuted data
Distances.pre_embedding = function(pbmc, pbmc_permuted, K, pre_embedding = 'pca') {
  distances <- distances::distances
  
  M = pbmc@reductions[[pre_embedding]]@cell.embeddings
  pre_embedding_distances = distances(M[,1:K])
  
  M_permuted = pbmc_permuted@reductions[[pre_embedding]]@cell.embeddings
  pre_embedding_distances_permuted = distances(M_permuted[,1:K])
  
  return(list("pre_embedding_distances"=pre_embedding_distances, "pre_embedding_distances_permuted"=pre_embedding_distances_permuted))
}

# Find distances under PCA and tSNE for both original and permuted matrices.
##default perplexity to match seurat
Distances.tSNE = function(pbmc,pbmc.permuted, K, perplexity_score = 40, pre_embedding = 'pca', check_duplicates = T, rerun = T)
{
  
  distances<-distances::distances
  if(rerun){
  pbmc = Seurat::RunTSNE(pbmc, seed.use=100, perplexity = perplexity_score, reduction = pre_embedding, do.fast = T, check_duplicates = check_duplicates)
  }
  tSNE_distances = distances(pbmc@reductions$tsne@cell.embeddings)

  pbmc.permuted <- Seurat::RunTSNE(pbmc.permuted, seed.use=100, perplexity = perplexity_score, reduction = pre_embedding, do.fast = T, check_duplicates = check_duplicates)
  tSNE_distances_permuted = distances(pbmc.permuted@reductions$tsne@cell.embeddings)
  
  results.PCA <- list("reduced_dim_distances"=tSNE_distances,"reduced_dim_distances_permuted"=tSNE_distances_permuted)
  
  return(results.PCA)
}

# Find distances under PCA and UMAP for both original and permuted matrices.
##default setting same as Seurat

Distances.UMAP = function(pbmc,pbmc.permuted, K, pre_embedding = 'pca', n = 30, m = 0.3, rerun = T) {
  distances <- distances::distances
  if(rerun){
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:K, seed.use = 100, reduction = pre_embedding, n.neighbors = n, min.dist = m)
}
  
  UMAP_distances = distances(pbmc@reductions$umap@cell.embeddings)
  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, dims = 1:K, seed.use = 100, reduction = pre_embedding, n.neighbors = n, min.dist = m)
  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)
  results.PCA   <- list("reduced_dim_distances" = UMAP_distances, "reduced_dim_distances_permuted" = UMAP_distances_permuted)
  
  return(results.PCA)
}


#calculate cell similarity scores 
Cell.Similarity = function(pre_embedding_distances,pre_embedding_distances_permuted,reduced_dim_distances,reduced_dim_distances_permuted,similarity_percent = .50) {
  numberselected = floor((dim(pre_embedding_distances)[2])*similarity_percent)
  rho_original  <- foreach::`%do%`(foreach::foreach(i=1:(dim(pre_embedding_distances)[2])), cor((reduced_dim_distances[i,order(pre_embedding_distances[i,])][2:(numberselected+1)]), (sort(reduced_dim_distances[i,])[2:(numberselected+1)])))
  
  rho_original = as.numeric(rho_original)
  # print("original done")
  
  rho_permuted <- foreach::`%do%`(foreach::foreach(i=1:(dim(pre_embedding_distances)[2])), cor((reduced_dim_distances_permuted[i, order(pre_embedding_distances_permuted[i,])][2:(numberselected+1)]), (sort(reduced_dim_distances_permuted[i,])[2:(numberselected+1)])))
  rho_permuted  = as.numeric(rho_permuted)
  # print("permuted done")
  
  similarity_results <- list("rho_original" = rho_original, "rho_permuted" = rho_permuted)
  return(similarity_results)
}

#classify cells based off null distribution of similarity scores
Cell.Classify = function(rho_original, rho_permuted, dubious_cutoff = 0.05, trustworthy_cutoff = 0.95) {
  y = seq(1,length(rho_original),length=length(rho_original))
  
  rho_trustworthy = quantile(rho_permuted,trustworthy_cutoff)
  rho_dubious     = quantile(rho_permuted,dubious_cutoff)
  
  dubious_cells     = which(rho_original < rho_dubious)
  trustworthy_cells = which(rho_original > rho_trustworthy)
  intermediate      = c(1:length(rho_original))
  intermediate      = intermediate[!intermediate%in% c(dubious_cells, trustworthy_cells)]
  
  ClassifiedCells.results <- list("dubious_cells" = dubious_cells,"trustworthy_cells"=trustworthy_cells,"intermediate_cells"=intermediate)
  
  return(ClassifiedCells.results)
}


#wrapper function
optimize = function(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                    n, m, perplexity, results.PCA, similarity_percent, dubious_cutoff,
                    trustworthy_cutoff, check_duplicates = T, rerun = T) {
  
  if(reduction.method == 'umap') {
    results <-Distances.UMAP(pbmc = input_data ,pbmc.permuted = input_data.permuted, K = K, pre_embedding = pre_embedding,
                              n = n, m = m, rerun = rerun)
  } else if (reduction.method == 'tsne') {
    results <- Distances.tSNE(pbmc = input_data ,pbmc.permuted = input_data.permuted, K = K, pre_embedding = pre_embedding,
                              perplexity_score = perplexity, check_duplicates = check_duplicates, rerun = rerun)
  }
  
  similarity_score <- Cell.Similarity(results.PCA$pre_embedding_distances, results.PCA$pre_embedding_distances_permuted,
                                      results$reduced_dim_distances, results$reduced_dim_distances_permuted, similarity_percent)
  ClassifiedCells  <- Cell.Classify(similarity_score$rho_original, similarity_score$rho_permuted, dubious_cutoff = dubious_cutoff,
                                    trustworthy_cutoff = trustworthy_cutoff)
  dub = ifelse(length(ClassifiedCells$dubious_cells)!=0, paste(ClassifiedCells$dubious_cells, sep = ',', collapse = ','), 'none')
  int = ifelse(length(ClassifiedCells$intermediate_cells)!=0, paste(ClassifiedCells$intermediate_cells, sep = ',', collapse = ','), 'none')
  trust = ifelse(length(ClassifiedCells$trustworthy_cells)!=0, paste(ClassifiedCells$trustworthy_cells, sep = ',', collapse = ','), 'none')
    
  return(c(length(ClassifiedCells$dubious_cells), dub, trust, int))
}


scDEED = function(input_data, K, n_neighbors = c(5, 20, 30, 40,50), min.dist = c(0.1, 0.4), similarity_percent = 0.5,reduction.method,
                  perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)), pre_embedding = 'pca', slot = 'scale.data', 
                  dubious_cutoff = 0.05, trustworthy_cutoff = 0.95, permuted = NA, check_duplicates = T, rerun = T, default_assay = 'active.assay'){
  

  if (is.na(permuted)) {
    print('Permuting data')
    input_data.permuted <- suppressMessages(Permuted(input_data, K = K,slot = slot, default_assay = default_assay))
    print('Permutation finished')
  } else {
    input_data.permuted = permuted
  }

  results.PCA = suppressMessages(Distances.pre_embedding(input_data, input_data.permuted, K=K, pre_embedding = pre_embedding))


  
  if (reduction.method == 'umap') {
    all_pairs <- expand.grid(n_neighbors, min.dist)
    if(nrow(all_pairs)>1){
      print('Estimating time for each hyperparameter setting...')
      start = Sys.time()
      original =  suppressMessages(foreach::`%do%`(foreach::foreach(n = all_pairs$Var1[1], m = all_pairs$Var2[1], .combine = "cbind"),
                                  optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                                           n = n, m = m, perplexity = NA, results.PCA = results.PCA,
                                           similarity_percent = similarity_percent,
                                           dubious_cutoff = dubious_cutoff,
                                           trustworthy_cutoff = trustworthy_cutoff, rerun = rerun)))
      end = Sys.time()
      time = end-start
      print('Estimated time per hyperparameter setting ')
      print(time)
      
      print('Estimated time of completion')
      print(Sys.time() + time*nrow(all_pairs))
      original = as.matrix(original)
    all_dub   <- suppressMessages(foreach::`%do%`(foreach::foreach(n = all_pairs$Var1[-1], m = all_pairs$Var2[-1], .combine = "cbind"),
                                 optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                                          n = n, m = m, perplexity = NA, results.PCA = results.PCA,
                                          similarity_percent = similarity_percent,
                                          dubious_cutoff = dubious_cutoff,
                                          trustworthy_cutoff = trustworthy_cutoff, rerun = rerun)))

    if(nrow(all_pairs)==2){
      all_dub=as.matrix(all_dub)
    } 
        all_dub = cbind(original, all_dub)
    } else{
      all_dub = suppressMessages(foreach::`%do%`(foreach::foreach(n = all_pairs$Var1, m = all_pairs$Var2, .combine = "cbind"),
                                       optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                                                n = n, m = m, perplexity = NA, results.PCA = results.PCA,
                                                similarity_percent = similarity_percent,
                                                dubious_cutoff = dubious_cutoff,
                                                trustworthy_cutoff = trustworthy_cutoff, rerun = rerun)))
      all_dub=as.matrix(all_dub)
      
    }
    all_dub = t(all_dub)
    colnames(all_dub) = c('number_dubious_cells', 'dubious_cells', 'trustworthy_cells', 'intermediate_cells')
    colnames(all_pairs) = c('n_neighbors', 'min.dist')
    dubious_number_UMAP <- cbind(all_pairs, all_dub)
    dub_para            <- data.frame(dubious_number_UMAP[, 1:3])
    dub_para_full =as.data.frame(dubious_number_UMAP)
    
   
  
  }
  

  
  if (reduction.method == "tsne") {
    
    if (is.null(input_data@reductions$tsne)) {
      input_data <- Seurat::RunTSNE(input_data, do.fast = T)

    }

if(length(perplexity)>1){
  print('Estimating time for each hyperparameter setting...')
  start = Sys.time()
  original  <-suppressMessages(foreach::`%do%`(foreach::foreach(p = perplexity[1], .combine = "cbind"),
                                                 optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                                                          n, m, perplexity = p, results.PCA = results.PCA,
                                                          similarity_percent = similarity_percent,
                                                          dubious_cutoff = dubious_cutoff, trustworthy_cutoff, check_duplicates = check_duplicates, rerun = rerun)))
  end = Sys.time()
  time = end-start
  print('Estimated time per hyperparameter setting ')
  print(time)
  original = as.matrix(original)
  
  print('Estimated time of completion')
  print(Sys.time() + time*length(perplexity))
    dubious_number_tSNE <-suppressMessages(foreach::`%do%`(foreach::foreach(p = perplexity[-1], .combine = "cbind"),
                                          optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                                                   n, m, perplexity = p, results.PCA = results.PCA,
                                                   similarity_percent = similarity_percent,
                                                   dubious_cutoff = dubious_cutoff, trustworthy_cutoff, check_duplicates = check_duplicates, rerun = rerun)))
    
    if(length(perplexity)==2){
      dubious_number_tSNE=as.matrix(dubious_number_tSNE)
    } 

  dubious_number_tSNE = cbind(original, dubious_number_tSNE)

      
} else{
  dubious_number_tSNE <-suppressMessages(foreach::`%do%`(foreach::foreach(p = perplexity, .combine = "cbind"),
                                                         optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                                                                  n, m, perplexity = p, results.PCA = results.PCA,
                                                                  similarity_percent = similarity_percent,
                                                                  dubious_cutoff = dubious_cutoff, trustworthy_cutoff, check_duplicates = check_duplicates, rerun = rerun)))
  dubious_number_tSNE = as.matrix(dubious_number_tSNE)
  }
  
    all_dub = t(dubious_number_tSNE)
    colnames(all_dub) = c('number_dubious_cells', 'dubious_cells', 'trustworthy_cells', 'intermediate_cells')

    dubious_number_tsne <- cbind(perplexity, all_dub)
    colnames(dubious_number_tsne)[1] = 'perplexity'
    dub_para            <- data.frame(dubious_number_tsne[, 1:2])
    if(length(dub_para)==1){
      dub_para = as.data.frame(t(dub_para))
    }
    dub_para_full =as.data.frame(dubious_number_tsne)
    
  }


  dub_para$number_dubious_cells = as.numeric(dub_para$number_dubious_cells)
  dub_para_full$number_dubious_cells = as.numeric(dub_para_full$number_dubious_cells)
  rownames(dub_para) = NULL
  rownames(dub_para_full) = NULL
  return(list(num_dubious = dub_para, full_results = dub_para_full))
}
