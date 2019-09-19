library(karyoploteR)
library(zoo)
library(animation)


set.seed(345666666)

#params
chr <- "chr1"
rand.min <- -1
rand.max <- 1

static.nimages <- 10


gg <- kp$chromosome.lengths

tt <- tileGenome(gg, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
jj <- cumsum(runif(length(tt)+2, min = rand.min, max=rand.max))
tt$y <- rollmean(jj, k=3)

other.data <- list()
for(i in seq_len(10)) {
  other.data[[i]] <- tileGenome(gg, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
  jj <- cumsum(runif(length(other.data[[i]])+2, min = rand.min, max=rand.max))
  other.data[[i]]$y <- rollmean(jj, k=3)
}

maxs <- unlist(lapply(other.data,function(x) max(abs(min(x$y)), abs(max(x$y)))))

#Prepare plots for animation
cols <- .karyoploter.colors$horizon$schemas$redblue6


max.abs <- ceiling(max(abs(min(tt$y)), abs(max(tt$y))))
ymax <- max.abs
ymin <- -1*max.abs

pp <- getDefaultPlotParams(1)
pp$bottommargin <- 30
pp$topmargin <- 30
pp$data1inmargin <- 10
pp$ideogramheight <- 15

chr.cex=2
axis.cex=1.6
axis.lwd=2
line.lwd=2

createPlots <- function() {
  
  # Img1
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
    #The DATA
    kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  }
  
  #Img2
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=1,num.breaks = 1, col=c("#FFFFFF00", "#AAAAAA"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0,num.breaks = 1, col=c("#AAAAAA", "#FFFFFF00"))
    
    #The DATA
    kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  }
  
  
  #Img3
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=1,num.breaks = 1, col=c("#FFFFFF00", cols[5]))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0,num.breaks = 1, col=c(cols[3], "#FFFFFF00"))
    
    #The DATA
    kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  }
  
  
  #Img4
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=1,num.breaks = 1, col=c("#FFFFFF00", cols[5]))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0,num.breaks = 1, col=c(cols[3], "#FFFFFF00"))
    
    #hlines
    kpAbline(kp, h=ymin+(ymax-ymin)/6*(1:5), col="#AAAAAA", ymin=ymin, ymax=ymax, lwd=2)
    
    #The DATA
    kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  }
  
  
  #Img5
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=1,num.breaks = 1, col=c("#FFFFFF00", cols[5]))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0,num.breaks = 1, col=c(cols[3], "#FFFFFF00"))
    
    #hlines
    kpAbline(kp, h=ymin+(ymax-ymin)/6*(1:6), col="white", ymin=ymin, ymax=ymax, lwd=2)
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
    #The DATA
    kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  }
  
  
  # #Img6
  # for(i in seq_len(static.nimages)) {
  #   kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
  #   
  #   kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=1,num.breaks = 1, col=c("#FFFFFF00", cols[5]))
  #   kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0,num.breaks = 1, col=c(cols[3], "#FFFFFF00"))
  #   
  #   #hlines
  #   kpAbline(kp, h=ymin+(ymax-ymin)/6*(1:6), col="white", ymin=ymin, ymax=ymax, lwd=2)
  #   
  #   kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
  #   kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
  #   
  #   #The DATA
  #   kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  # }
  # 
  
  #Img7
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    
    #Pos
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[5], "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[6], "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[7]))
    
    #Neg
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5-0.5/3*(1-1), r1=0.5-0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", cols[3], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5-0.5/3*(2-1), r1=0.5-0.5/3*2, col=c("#FFFFFF00", cols[2], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5-0.5/3*(3-1), r1=0.5-0.5/3*3, col=c(cols[1], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    
    #hlines
    kpAbline(kp, h=ymin+(ymax-ymin)/6*(1:6), col="white", ymin=ymin, ymax=ymax, lwd=2)
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
    #The DATA
    kpLines(kp, tt, ymin=ymin, ymax=ymax, r0=0, r1=1, lwd=line.lwd)
  }
  
  #Img8
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    
    #Pos
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[5], "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[6], "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[7]))
    
    #Neg
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5-0.5/3*(1-1), r1=0.5-0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", cols[3], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5-0.5/3*(2-1), r1=0.5-0.5/3*2, col=c("#FFFFFF00", cols[2], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5-0.5/3*(3-1), r1=0.5-0.5/3*3, col=c(cols[1], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
  }
  
  
  
  
  
  
  #Img9-13
  nimages <- 10 
  
  for(i in seq_len(nimages+1)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    #Pos
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[5], "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[6], "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[7]))
    
    #Neg
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+(-0.5+1/nimages*(i-1))/3*(1-1), r1=0.5+(-0.5+1/nimages*(i-1))/3*1, col=c("#FFFFFF00", "#FFFFFF00", cols[3], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+(-0.5+1/nimages*(i-1))/3*(2-1), r1=0.5+(-0.5+1/nimages*(i-1))/3*2, col=c("#FFFFFF00", cols[2], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+(-0.5+1/nimages*(i-1))/3*(3-1), r1=0.5+(-0.5+1/nimages*(i-1))/3*3, col=c(cols[1], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
    
  }
  
  #Img13
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    #Pos
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[5], "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[6], "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[7]))
    
    #Neg
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", cols[3], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", cols[2], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c(cols[1], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    
    kpAxis(kp, ymin=ymin, ymax=ymax, cex=axis.cex, lwd=axis.lwd)
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
  }
  
  
  #Img14
  for(i in seq_len(static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    #Pos
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[5], "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[6], "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[7]))
    
    #Neg
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", cols[3], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1), r1=0.5+0.5/3*2, col=c("#FFFFFF00", cols[2], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1), r1=0.5+0.5/3*3, col=c(cols[1], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
  }
  
  
  
  
  #Img15-20
  nimages <- 9
  
  for(i in seq_len(nimages+1)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    
    #Pos
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+(0.5/3*1), col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[5], "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1)-(0.5/3*1)/nimages*(i-1), r1=0.5+0.5/3*2-(0.5/3*1)/nimages*(i-1), col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[6], "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1)-(0.5/3*2)/nimages*(i-1), r1=0.5+0.5/3*3-(0.5/3*2)/nimages*(i-1), col=c("#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", cols[7]))
    
    #Neg
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(1-1), r1=0.5+0.5/3*1, col=c("#FFFFFF00", "#FFFFFF00", cols[3], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(2-1)-(0.5/3*1)/nimages*(i-1), r1=0.5+0.5/3*2-(0.5/3*1)/nimages*(i-1), col=c("#FFFFFF00", cols[2], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5+0.5/3*(3-1)-(0.5/3*2)/nimages*(i-1), r1=0.5+0.5/3*3-(0.5/3*2)/nimages*(i-1), col=c(cols[1], "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"))
    
    kpAbline(kp, h=0, ymin=ymin, ymax=ymax)
  }
  
  
  
  #Img21
  for(i in seq_len(2*static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
      
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
  }
  
  #Img22-25
    margin <- 0.01
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
    kpPlotHorizon(kp, other.data[[1]], ymin=-1*maxs[1], ymax=maxs[1], r0=0.5-0.5/3-1*margin, r1=0.5-1*margin, col=cols)
    
   
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
    kpPlotHorizon(kp, other.data[[1]], ymin=-1*maxs[1], ymax=maxs[1], r0=0.5-0.5/3-1*margin, r1=0.5-1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[8]], ymin=-1*maxs[8], ymax=maxs[8], r0=0.5+0.5/3+1*margin, r1=0.5+0.5/3*2+1*margin, col=cols)

    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
    kpPlotHorizon(kp, other.data[[1]], ymin=-1*maxs[1], ymax=maxs[1], r0=0.5-0.5/3-1*margin, r1=0.5-1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[8]], ymin=-1*maxs[8], ymax=maxs[8], r0=0.5+0.5/3+1*margin, r1=0.5+0.5/3*2+1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[3]], ymin=-1*maxs[3], ymax=maxs[3], r0=0.5-0.5/3*2-2*margin, r1=0.5-0.5/3-2*margin, col=cols)

    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
    kpPlotHorizon(kp, other.data[[1]], ymin=-1*maxs[1], ymax=maxs[1], r0=0.5-0.5/3-1*margin, r1=0.5-1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[8]], ymin=-1*maxs[8], ymax=maxs[8], r0=0.5+0.5/3+1*margin, r1=0.5+0.5/3*2+1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[3]], ymin=-1*maxs[3], ymax=maxs[3], r0=0.5-0.5/3*2-2*margin, r1=0.5-0.5/3-2*margin, col=cols)
    kpPlotHorizon(kp, other.data[[9]], ymin=-1*maxs[9], ymax=maxs[9], r0=0.5+0.5/3*2+2*margin, r1=0.5+0.5/3*3+2*margin, col=cols)

    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
    kpPlotHorizon(kp, other.data[[1]], ymin=-1*maxs[1], ymax=maxs[1], r0=0.5-0.5/3-1*margin, r1=0.5-1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[8]], ymin=-1*maxs[8], ymax=maxs[8], r0=0.5+0.5/3+1*margin, r1=0.5+0.5/3*2+1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[3]], ymin=-1*maxs[3], ymax=maxs[3], r0=0.5-0.5/3*2-2*margin, r1=0.5-0.5/3-2*margin, col=cols)
    kpPlotHorizon(kp, other.data[[9]], ymin=-1*maxs[9], ymax=maxs[9], r0=0.5+0.5/3*2+2*margin, r1=0.5+0.5/3*3+2*margin, col=cols)
    kpPlotHorizon(kp, other.data[[7]], ymin=-1*maxs[7], ymax=maxs[7], r0=0.5-0.5/3*3-3*margin, r1=0.5-0.5/3*2-3*margin, col=cols)
    
  #img Final
  for(i in seq_len(4*static.nimages)) {
    kp <- plotKaryotype(chr=chr, plot.params = pp, cex=chr.cex)
    kpPlotHorizon(kp, tt, ymin=ymin, ymax=ymax, r0=0.5, r1=0.5+(0.5/3), col=cols)
    kpPlotHorizon(kp, other.data[[1]], ymin=-1*maxs[1], ymax=maxs[1], r0=0.5-0.5/3-1*margin, r1=0.5-1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[8]], ymin=-1*maxs[8], ymax=maxs[8], r0=0.5+0.5/3+1*margin, r1=0.5+0.5/3*2+1*margin, col=cols)
    kpPlotHorizon(kp, other.data[[3]], ymin=-1*maxs[3], ymax=maxs[3], r0=0.5-0.5/3*2-2*margin, r1=0.5-0.5/3-2*margin, col=cols)
    kpPlotHorizon(kp, other.data[[9]], ymin=-1*maxs[9], ymax=maxs[9], r0=0.5+0.5/3*2+2*margin, r1=0.5+0.5/3*3+2*margin, col=cols)
    kpPlotHorizon(kp, other.data[[7]], ymin=-1*maxs[7], ymax=maxs[7], r0=0.5-0.5/3*3-3*margin, r1=0.5-0.5/3*2-3*margin, col=cols)
  }
  
}


saveGIF(createPlots(), movie.name = "horizon.animation.gif", interval=0.1, ani.width=1000, ani.height=800)



