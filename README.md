# Examples of plots created with karyoploteR

This repo contains a number of small examples on how to use the Bioconductor package [karyoploteR](http://bioconductor.org/packages/karyoploteR) 
to create karyoplots, that is, the chromosome ideograms accompained by any arbitrary
data.

**IMPORTANT NOTE:** At the moment karyoploteR is only available in the _devel_ version 
of Bioconductor. Therefore it will only be available on R-devel instances. Please refer to the 
documentation on [how to use Bioconductor devel](https://www.bioconductor.org/developers/how-to/useDevel/) 
for more information.

The first part is a tutorial-like series of simple examples showing the different parts
of karyoploteR and how to use them. The second part contains more involved and complete
stand-alone examples. 

In all cases, each example is self contained and it should be possible to directly copy 
and run the R code.

## Tutorial

* [Example 1 - Plot ideograms](Examples/Tutorial/CreateIdeogram/CreateIdeogram.md)

    Create empty ideograms (with no data plotted) for different organisms
    
     
* [Example 2 - Filter and reorder chromosomes](Examples/Tutorial/FilterChromosomes/FilterChromosomes.md)

    Create ideograms of a subset of chromosomes and plot them in any order
    
* [Example 3 - Using custom genomes](Examples/Tutorial/CustomGenomes/CustomGenomes.md)
    
    Create ideograms using your own custom genomes including specifying your own cytobands

* [Example 4 - Adding base numbers and cytoband labels](Examples/Tutorial/BaseNumbersAndBandNames/BaseNumbersAndBandNames.md)
    
    Add a base numbering guide to the ideograms and label the cytobands with their names



* More to come soon...    
    
## Complete Examples

* [Example 1 - Multiple data types](Examples/CompleteExamples/MultipleDataTypes/MultipleDataTypes.md)

    A plot with multiple data types plotted together. 
    
    ![Multiple Data Types Figure](Examples/CompleteExamples/MultipleDataTypes/figure/Figure-1.png?raw=true "Multiple Data Types")


*** 

For more extensive documentation on keryoploteR, please refer to the [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf) and the [reference manual](http://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf)