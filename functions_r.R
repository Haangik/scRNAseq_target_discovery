get_tumor_normal_expr<-function(gene, margin_data){
  res<-margin_data %>% dplyr::filter(Gene==gene) %>%
    dplyr::select(Gene, mean_expr_cancer, mean_expr_normal, log2FC_cancer_vs_normal) %>%
    dplyr::slice_max(mean_expr_cancer, n = 1, with_ties = FALSE)
  
  if(nrow(res)==0){
    res<-data.frame(Gene=gene, mean_expr_cancer=NA, mean_expr_normal=NA, log2FC_cancer_vs_normal=NA)
    return(res)
  }
  return(res)
}

compute_f1_score<-function(adata, g1, g2){
  res<-py$compute_positivity(adata=adata, gene1=g1, gene2=g2)[c(4,5)]
  score=2*prod(res)/sum(res)
  return(score)
}
