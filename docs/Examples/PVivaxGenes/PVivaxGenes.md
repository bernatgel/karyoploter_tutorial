---
layout: example
label: PVivaxGenes
title: PVivax Genes
---




## Plot Genes in Plasmodium Vivax Genome

In this example we'll get the genome of Plasmodium Vivax and it's genes from 
a _gff_ file downloaded from [PlasmoDB](http://plasmodb.org). In this case,
the _gff_ file contains the lengths of the chromosomes, but it could be
extracted from any other source.

Since we want to use the chromosome lengths contained in the gff file header, 
we'll use `readLines` to read part of the header (30 lines), where we 
know that the genome information is located (you have to take a look at the file
first to know this and take into account that most gff's do not have this 
information). After that, we'll extract the information by splitting the lines
and selecting the appropiate parts. This is not the most elegant solution to 
this problem, but it works.



```r
library(karyoploteR)

gff.file <- "http://plasmodb.org/common/downloads/Current_Release/PvivaxP01/gff/data/PlasmoDB-33_PvivaxP01.gff"

header.lines <- readLines(gff.file, n = 30)
```

```
## Error in file(con, "r"): cannot open the connection to 'http://plasmodb.org/common/downloads/Current_Release/PvivaxP01/gff/data/PlasmoDB-33_PvivaxP01.gff'
```

```r
#The lines with the standard chromosomes start with "##sequence-region PvP01".
#Select them.
ll <- header.lines[grepl(header.lines, pattern = "##sequence-region PvP01")]
```

```
## Error in eval(expr, envir, enclos): object 'header.lines' not found
```

```r
#split them by space, and create a data.frame
gg <- data.frame(do.call(rbind, strsplit(ll, split = " ")))
```

```
## Error in strsplit(ll, split = " "): object 'll' not found
```

```r
gg[,3] <- as.numeric(as.character(gg[,3]))
```

```
## Error in eval(expr, envir, enclos): object 'gg' not found
```

```r
gg[,4] <- as.numeric(as.character(gg[,4]))
```

```
## Error in eval(expr, envir, enclos): object 'gg' not found
```

```r
#and create a GRanges with the information
PvP01.genome <- toGRanges(gg[,c(2,3,4)])
```

```
## Error in is(A, "GRanges"): object 'gg' not found
```

```r
PvP01.genome
```

```
## Error in eval(expr, envir, enclos): object 'PvP01.genome' not found
```

With this, we can already plot the genome of Plasmodium Vivax using the function
`plotKaryotype`.


```r
kp <- plotKaryotype(genome=PvP01.genome)
```

```
## Error in is(genome, "GRanges"): object 'PvP01.genome' not found
```

Since we used a custom genome and provided no cytobands information, 
`plotKaryotype` created a representation of the genome with gray cromosomes. 
Actually, we could substitute the gray boxe by gray lines to simplify them a bit
calling `kpAddCytobandsAsLine` instead of the default `kpAddCytobands` by setting
`ideogram.plotter=NULL`.



```r
kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL)
```

```
## Error in is(genome, "GRanges"): object 'PvP01.genome' not found
```

```r
kpAddCytobandsAsLine(kp)
```

```
## Error in kpAddCytobandsAsLine(kp): object 'kp' not found
```

Now, to plot the genes on the genome, we need to download them. We'll use the 
same _gff_ file, but now we'll use the `import` function from 
[rtracklayer](http://bioconductor.org/packages/rtracklayer/) to download and 
read the file content into a GRanges object. The function is automatically 
imported when loading karyoploteR. 

Since the gff file contains featues of many different typesm we'll filter it and
keep only the genes, since this is what we want to plot.



```r
features <- import(gff.file)
```

```
## Error in download.file(resource(con), destfile): cannot open URL 'http://plasmodb.org/common/downloads/Current_Release/PvivaxP01/gff/data/PlasmoDB-33_PvivaxP01.gff'
```

```r
table(features$type)
```

```
## Error in features$type: object of type 'closure' is not subsettable
```

```r
genes <- features[features$type=="gene"]
```

```
## Error in features$type: object of type 'closure' is not subsettable
```

With that, we can plot the genes in the chromosomes. In this case, since we'll 
plot more than 6800 genes, we'll plot them as regions, with no indication of the
gene name or id.


```r
kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL)
```

```
## Error in is(genome, "GRanges"): object 'PvP01.genome' not found
```

```r
kpAddCytobandsAsLine(kp)
```

```
## Error in kpAddCytobandsAsLine(kp): object 'kp' not found
```

```r
kpPlotRegions(kp, data=genes)
```

```
## Error in methods::is(karyoplot, "KaryoPlot"): object 'kp' not found
```

We can not see a lot in the plot, but we can see that some of the genes are 
above the rest, because they overlap another gene. Since in this case we'd prefer
all genes to be in a single line, we can add `avoid.overlapping = FALSE` to
the call to `kpPlotRegions`. In addition, we can add a second data panel below
the choromosomes (using `plot.type=2`) and plot the forward strand genes above 
the ideogram and the reverse strand genes below the ideogram. To do that we'll 
need to filter the genes by _strand_ and explicitly ask for the genes in the 
negative strand to be plotted in the second data panel (below the ideogram) 
with `data.panel=2`.


```r
kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL, plot.type=2)
```

```
## Error in is(genome, "GRanges"): object 'PvP01.genome' not found
```

```r
kpAddCytobandsAsLine(kp)
```

```
## Error in kpAddCytobandsAsLine(kp): object 'kp' not found
```

```r
kpPlotRegions(kp, data=genes[strand(genes)=="+"], avoid.overlapping = FALSE)
```

```
## Error in methods::is(karyoplot, "KaryoPlot"): object 'kp' not found
```

```r
kpPlotRegions(kp, data=genes[strand(genes)=="-"], avoid.overlapping = FALSE, data.panel=2)
```

```
## Error in methods::is(karyoplot, "KaryoPlot"): object 'kp' not found
```

With that we can see some more interesting things, such as stetches of the
genome with only forward or only reverse genes, but the plot is still difficult 
to read. We'd like to separate the chromosomes between them and to do that, we 
know we'll have to modify the `plot.params`, but what plot param exactly? 
We can plot them all and see what we want to change.


```r
plotDefaultPlotParams(plot.type=2)
```

![plot of chunk Figure5](images//Figure5-1.png)

Looking at the figure, we can see that the parameters affecting the distance 
between chromosomes are _data1outmargin_ and _data2outmargin_ and that they have
a default value of 20. We'll increase them to 100 to get a good separation 
between chromosomes and we'll also add a main title (and some more margin for it)
and colors to differentiate the plus and minus genes.


```r
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
kp <- plotKaryotype(genome=PvP01.genome, ideogram.plotter = NULL, plot.type=2, plot.params = pp)
```

```
## Error in is(genome, "GRanges"): object 'PvP01.genome' not found
```

```r
kpAddCytobandsAsLine(kp)
```

```
## Error in kpAddCytobandsAsLine(kp): object 'kp' not found
```

```r
kpAddMainTitle(kp, "Plasmodium Vivax - PvP01 with genes", cex=2)
```

```
## Error in kpAddMainTitle(kp, "Plasmodium Vivax - PvP01 with genes", cex = 2): object 'kp' not found
```

```r
kpPlotRegions(kp, data=genes[strand(genes)=="+"], avoid.overlapping = FALSE, col="deepskyblue")
```

```
## Error in methods::is(karyoplot, "KaryoPlot"): object 'kp' not found
```

```r
kpPlotRegions(kp, data=genes[strand(genes)=="-"], avoid.overlapping = FALSE, col="gold", data.panel=2)
```

```
## Error in methods::is(karyoplot, "KaryoPlot"): object 'kp' not found
```

```r
kpAddLabels(kp, "strand +", cex=0.8, col="#888888")
```

```
## Error in kpAddLabels(kp, "strand +", cex = 0.8, col = "#888888"): object 'kp' not found
```

```r
kpAddLabels(kp, "strand -", data.panel=2, cex=0.8, col="#888888")
```

```
## Error in kpAddLabels(kp, "strand -", data.panel = 2, cex = 0.8, col = "#888888"): object 'kp' not found
```
