# defining global variables and functions to appease R CMD Check

utils::globalVariables(
  names = c(
    ".",".data", "baseMean","biotype","chr_start","chromosome","col_sum","colname","count", "detection_threshold","fraction_of_counts",
    "fraction_of_genes","gc_content","gene_biotype","gene_id","gene_mean","gene_name","gene_sd","log_count", "log2FoldChange",
    "number_of_genes","padj","PC_rank","pvalue", "rowData","rowname","sample_id","seqnames","significant","start","svalue",
    "total_count","tx_id","var_exp","x","y","ylim"
  )
)
