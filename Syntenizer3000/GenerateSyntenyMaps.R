if(!require(MASS)){
  install.packages("MASS", repos = "http://cloud.r-project.org/")
  library(MASS)
}
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1)
{
  stop("When running this script please provide path of folder containing map data (the output folder of Syntenizer3000) as an argument.", call.=FALSE)
} else
{
  setwd(args[1])
  if (file.exists("gene_groups.csv"))
  {
    groups<-read.csv(file="gene_groups.csv", sep=":", header = F)[,1]
    print("Rendering synteny maps, please wait...")
    progress <- -1
    for (group in 1:length(groups)){
      filename <- paste(groups[group], "_synteny_chart.csv", sep="")
      map = read.csv(file = filename, dec=",", sep=";", header=T, check.names = F)
      strains = colnames(map)
      
      pdf(paste(groups[group], "_synteny_chart.pdf", sep=""), width=40, height=20)
      
      par(mar=c(10,7,0,0))
      par(las=2, cex=1.2, lwd=1.0)
      parcoord(map, lty = 1, main="", col= c("#00000000", "#00000000", "black", sample(rainbow(nrow(map)-3))), lwd=3)
      axis(2, at=seq(0,1, 1/42), labels=c("Upstream", seq(-20,-1,1 ), "GENE GROUP", paste(paste0("+", as.character(seq(1,20,1 )))) ,"Downstream"))

      #par(mar=c(10,7,0,0))
      #par(las=2, cex=0.5, lwd=0.25)
      #parcoord(map, lty = 1, main="", col= c("#00000000", "#00000000", "black", sample(rainbow(nrow(map)-3))))
      #axis(2, at=seq(0,1, 1/42), labels=c("Upstream", seq(-20,-1,1 ), "GENE GROUP", paste(paste0("+", as.character(seq(1,20,1 )))) ,"Downstream"))
      
      dev.off()
      file.remove(filename)
      
      if (progress != as.integer( (group*10)/length(groups)))
      {
        progress = as.integer((group*10)/length(groups))
        print(paste(progress*10, "%", sep="") )
      }
    }
    print("Rendering complete.")
  }
  else
  {
    stop("Could not locate gene_groups.csv", call.=FALSE)
  }
}

