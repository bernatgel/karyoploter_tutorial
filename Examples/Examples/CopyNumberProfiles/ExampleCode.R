library(RTCGA)
library(karyoploteR)


(cohorts <- infoTCGA() %>% 
  rownames() %>% 
  sub("-counts", "", x=.))

cohort <- "BLCA"
releaseDate <- "2015-11-01"

downloadTCGA(cancerTypes = cohort, destDir = tempdir(), date = releaseDate, dataSet = "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3" )

cnv.file <- list.files(tempdir(), pattern = "cnv", recursive = TRUE, full.names = TRUE)
cnv.data <- read.table(cnv.file, sep="\t", header=TRUE, stringsAsFactors = FALSE)

hist(cnv.data$Segment_Mean)

cnv.data <- setNames(cnv.data[, c(2,3,4,6,1,5)], c("chr", "start", "end", "y", "sample","probes"))
cnv.data <- toGRanges(cnv.data)
seqlevelsStyle(cnv.data) <- "UCSC"
cnv.data <- split(cnv.data, cnv.data$sample)
unlink(cnv.file)


kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL)
kpAddCytobandsAsLine(kp)
len <- length(cnv.data)
for(i in seq_along(cnv.data)) {
  kpHeatmap(kp, data = cnv.data[[i]], colors = c('blue','white', "white", 'red'), r0=(i-1)*1/len, r1=i*1/len)
}
#kpPlotMarkers(kp, chr = "chr9", x=21967752, labels="CDK2NA")
