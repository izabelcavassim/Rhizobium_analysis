if(!require(circlize)){
  install.packages("circlize", repos = "http://cloud.r-project.org/")
  library(circlize)
}
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1)
{
  stop("When running this script please provide path of folder containing chart data (the output folder of Syntenizer3000) as an argument.", call.=FALSE)
} else
{
  setwd(args[1])
  if (file.exists('strains.csv'))
  {
    strains = read.csv(file = 'strains.csv', sep=";")
    print("Rendering strain charts, please wait...")
    progress <- -1
    for(s in 1:nrow(strains)){
      if (progress != as.integer( (s*10)/nrow(strains)))
      {
        progress = as.integer((s*10)/nrow(strains))
        print(paste(progress*10, "%", sep="") )
      }
      filename_prefix <- strains$ID[s]
      pdf(paste(filename_prefix,'_chart.pdf',sep=""), width=50, height=50)
      
      filename.contigs <- paste(filename_prefix,'_contigs.csv',sep="")
      filename.lanes <- paste(filename_prefix,'_lanes.csv',sep="")
      filename.paralogs <- paste(filename_prefix,'_paralogs.csv',sep="")
      contigs = read.csv(file = filename.contigs, dec=".", sep=";", stringsAsFactors = F, colClasses = c("ID"="character"))
      lanes = read.csv(file = filename.lanes, dec=".", sep=";", stringsAsFactors = F, colClasses = c("ID"="character"))
      paralogs = read.csv(file = filename.paralogs, dec=".", sep=";", stringsAsFactors = F, colClasses = c("ID1"="character", "ID2"="character"))

      lanes <- lanes[order(lanes$ID, lanes$Location),]

      circos.clear()
      circos.par(cell.padding = c(0.0, 0, 0.0, 0), track.height = 0.25, start.degree = 90, gap.degree = 0.5, track.margin=c(0,0.01))
      circos.initialize(factors = unique(contigs$ID), xlim = t(matrix(data = contigs$Bound, ncol = nrow(contigs)/2, nrow = 2)))
      circos.trackPlotRegion(factors = as.vector.factor(unique(contigs$ID)), ylim=c(0,1), panel.fun = function(x, y) {
        theta <- circlize(CELL_META$xcenter, 1.3)[1, 1] %% 360
        dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
        aa <- c(1.15,0.5)
        if (theta < 90 || theta > 270) aa <- c(-0.15,0.5)
        circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2],
                    CELL_META$sector.index, facing = dd, niceFacing = FALSE
                    , adj = aa
        )
        if(CELL_META$xlim[2] > 100000)
          circos.axis(labels=paste(seq(from=1,to=CELL_META$xlim[2]/100000,by=1), "00 Kb", sep=""), labels.cex=0.6, direction = "outside", major.at=seq(from=100000,to=CELL_META$xlim[2],by=100000), 
                      minor.ticks=1, labels.away.percentage = 0.15)
        else if (CELL_META$xlim[2] > 5000)
          circos.axis(labels=paste(round(CELL_META$xlim[2]/1000, 1), "Kb", sep=" "), labels.cex=0.6, direction = "outside", major.at=CELL_META$xlim[2], 
                      minor.ticks=0, labels.away.percentage = 0.15)
      })
      
      contig.colours<-contigs[ which(contigs$Bound!=0), ]
      for(i in 1:nrow(contig.colours)){
        contig <- as.data.frame(lanes[ which(lanes$ID==contig.colours$ID[i]), ])
        circos.trackLines(contig$ID, contig$Location, contig$Synteny, col = as.character( contig.colours$Colour[i]), border = "black", lwd=0.25, area=TRUE)
      }
      
      gene.colours <- unique(lanes$Colour)
      for (colour in gene.colours)
        if (!is.na(colour))
        {
          lines <- lanes[which(lanes$Colour == colour), ]
          circos.trackLines(lines$ID, lines$Location, lines$Synteny, col = colour, type = "h", lwd=1.5)
        }
      
      ## GC3s plot
      lanes.gc3 <- lanes[which(!is.na(lanes$GC3s)),] 
      if (nrow(lanes.gc3) > 0)
      {
        circos.par(cell.padding = c(0.0, 0, 0.0, 0), track.height = 0.20)
        circos.trackPlotRegion(factors = as.vector.factor(unique(contigs$ID)), ylim=c(0,1))
        
        circos.trackLines(contigs$ID, contigs$Bound, rep(0.5, nrow(contigs)), col = "grey")
        circos.trackLines(lanes$ID, lanes$Location, lanes$GC3s, col = "black", lwd=1, type="h", baseline=mean(lanes$GC3s))
        
        for (colour in gene.colours)
          if (!is.na(colour))
          {
            lines <- lanes[which(lanes$Colour == colour), ]
            circos.trackLines(lines$ID, lines$Location, lines$GC3s, col = colour, type = "h", lwd=1.5)
          }
      }
      
      ## Abundance plot
      circos.par(cell.padding = c(0.0, 0, 0.0, 0), track.height = 0.10)
      circos.trackPlotRegion(factors = as.vector.factor(unique(contigs$ID)), ylim=c(0,1))
      
      for(i in 1:nrow(contig.colours))
        circos.updatePlotRegion(sector.index = contig.colours$ID[i], bg.col = "grey")
      
      circos.trackLines(lanes$ID, lanes$Location, lanes$Abundance, col = "black", border = "black", lwd=0.25, area=TRUE)
      
      for (colour in gene.colours)
        if (!is.na(colour))
        {
          lines <- lanes[which(lanes$Colour == colour), ]
          circos.trackLines(lines$ID, lines$Location, lines$Abundance, col = colour, type = "h", lwd=1.5)
        }
      ## Unique genes plot
      circos.par(cell.padding = c(0.0, 0, 0.0, 0), track.height = 0.015, track.margin=c(0,0))
      circos.trackPlotRegion(factors = as.vector.factor(unique(contigs$ID)), ylim=c(0,1))
      for(i in 1:nrow(contig.colours))
        circos.updatePlotRegion(sector.index = contig.colours$ID[i], bg.col = "grey")
      
      uniquegenes <- lanes[ which(lanes$Abundance == 0), ]
      if (nrow(uniquegenes) > 0)
        circos.trackLines(uniquegenes$ID, uniquegenes$Location, matrix(data = 1, nrow = nrow(uniquegenes), ncol = 1),  col = "red", type = "h", lwd=0.5)
  
      #Paralogs plot
      if (nrow(paralogs) > 0)
        for(i in 1:nrow(paralogs))
          circos.link(paralogs$ID1[i], paralogs$Location1[i], as.character(paralogs$ID2[i]), paralogs$Location2[i], col = rand_color(1), lwd=0.5)
      
      dev.off()
      file.remove(filename.contigs)
      file.remove(filename.lanes)
      file.remove(filename.paralogs)
    }
    file.remove('strains.csv')
    print("Rendering complete.")
  } else
  {
    stop("Could not locate strains.csv", call.=FALSE)
  }
}