args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1)
{
  stop("When running this script please provide path of folder containing chart data (the output folder of Syntenizer3000) as an argument.", call.=FALSE)
} else
{
  setwd(args[1])
  if (file.exists('gene_relation_matrix.csv'))
  {
    grm <- as.matrix(read.csv('gene_relation_matrix.csv', sep = ';', header = T, stringsAsFactors=FALSE, check.names = FALSE)[,-1])
    rownames(grm)<-colnames(grm)
    
    pdf("Heatmap.pdf", width=15, height=15)
    heatmap(grm, symm = TRUE)
    dev.off()

    print("Rendering Heatmap.pdf complete.")
    
  } else
  {
    stop("Could not locate gene_relation_matrix.csv", call.=FALSE)
  }
}
