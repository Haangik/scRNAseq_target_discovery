# Example code
rm(list=ls())
library(dplyr)
library(reticulate)
library(ggplot2)
library(pbapply)
library(readr)

# Tutorial for Lung EGFR Post-TKI patient data
proj_name="Lung_primary_post-EGFR-TKI"

# if there is conda environment, activate the environment
#use_condaenv("haangik_py311")
py_run_file("LoadData.py")

# loading R functions
source("functions_r.R")


# Desired genes
gene_cand<-c("CEACAM5", "CEACAM6", "CLDN3", "LGR5")

# Normal-to-tumor margin analysis
result_normal_tumor_margin=py$compute_expression_margins(adata=py$adata_lung_primary_all,
                                     celltype_col = "cell_type_reannotation",
                                     gene_symbol_col="GeneSymbol",
                                     batch_size=as.integer(3000),
                                     immune_labels=c("T cell", "NK", "Plasma cell", "B cell",
                                                     "Mast", "Endothelial cell", "Myeloid cell")) %>% filter(!is.na(Gene))

# Positivity analysis
result_positivity=py$summarize_positivity(adata_tumor=py$adata_lung_primary_tumor,
                                          adata_normal=py$adata_lung_normal,
                                          genes=gene_cand,
                                          target1="CEACAM5")
## Calculating F1 score of two target genes
result_positivity$F1_Tumor<-pbsapply(gene_cand1, FUN=function(g2){
  score<-tryCatch(
    compute_f1_score(adata=py$adata_lung_primary_tumor, g1="CEACAM5", g2=g2),
    error=function(e) NA_real_
  )
  return(score)
})
result_positivity$F1_Normal<-pbsapply(gene_cand1, FUN=function(g2){
  score<-tryCatch(
    compute_f1_score(adata=py$adata_lung_normal, g1="CEACAM5", g2=g2),
    error=function(e) NA_real_
  )
  return(score)
})


# Dot plots of the desired genes
dotplot_dir=".../dotplot_dir/"
setwd(dotplot_dir)
dir.create(proj_name)
setwd(proj_name)
cell_types <- unique(py$adata$obs[["cell_type_reannotation"]])

# plotting dot plots
pbsapply(gene_cand, FUN=function(g){
  py$symbol=g
  py_run_string("
dp = sc.pl.dotplot(
    adata,
    var_names=symbol,
    groupby='cell_type_reannotation',
    standard_scale='var',
    color_map='Reds',
    dot_max=1.0,
    show=False
)

plt.savefig(f'{symbol}.png', dpi=300, bbox_inches='tight')
plt.close()
gc.collect()
")
return()
})

# Violin plots of the desired genes
violinplot_dir=".../violinplot_dir"
setwd(violinplot_dir)
dir.create(proj_name)
setwd(proj_name)

pbsapply(gene_cand, FUN=function(g){
  gene=g
  cell_types<-factor(c("Cancer cell", "Epithelial cell", gene_cell_stats$immune_celltype_max[i]),
                     levels=c("Cancer cell", "Epithelial cell", gene_cell_stats$immune_celltype_max[i]))
  df_expression<-do.call(rbind, lapply(cell_types, FUN=function(ct){
    vec <- py$cell_specific_gene_expression(
      adata = py$adata,
      gene = gene,
      celltype_col = "cell_type_reannotation",
      celltype_label = ct,
      gene_symbol_col = "GeneSymbol"
    )
    res<-data.frame(
      expression = as.numeric(vec),
      cell_type = ct
    )
    return(res)
  }))
  
  
  max_expr <- max(df_expression$expression, na.rm = TRUE)
  label_y <- max_expr * 1.05
  summary_df <- df_expression %>%
    group_by(cell_type) %>%
    summarise(
      mean_expr = mean(expression, na.rm = TRUE),
      median_expr = median(expression, na.rm = TRUE)
    ) %>%
    mutate(label_y = label_y) 
  
  p=ggplot(df_expression, aes(x = cell_type, y = expression, fill = cell_type)) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      color = NA,
      alpha = 0.5,
      width = 1.0
    ) +
    geom_boxplot(
      width = 0.05,
      outlier.size = 0.5,
      color = "black",
      fill = "white"
    ) +
    geom_text(
      data = summary_df,
      aes(
        x = cell_type,
        y = label_y,
        label = paste0("Mean: ", round(mean_expr, 2), "\nMedian: ", round(median_expr, 2))
      ),
      inherit.aes = FALSE,
      size = 3.5,
      vjust = 0
    ) +
    scale_y_sqrt(
      limits = c(0, label_y * 1.05)  # y축 범위 고정
    ) +
    theme_bw(base_size = 14) +
    labs(
      x = NULL,
      y = paste0("Normalized Expression (CP10K)"),
      title = paste0(gene, " expression")
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "none"
    )
  filename=paste0(gene, ".png", collapse = "")
  suppressWarnings(ggsave(filename, plot=p, width=8, height=6.5))
  py_run_string("gc.collect()")
  return()
})

