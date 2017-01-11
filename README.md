# Examples of plots created with karyoploteR

This repo contains a number of small examples on how to use the Bioconductor package [karyoploteR](http://bioconductor.org/packages/release/bioc/html/karyoploteR.html) 
to create karyoplots, that is, the chromosome ideograms accompained by any arbitrary
data.

The first part is a tutorial-like series of simple examples showing the different parts
of karyoploteR and how to use them. The second part contains more involved and complete
stand-alone examples. 

In all cases, each example is self contained and it should be possible to directly copy 
and run the R code.

## Tutorial

* [Example 1 - Create Random Regions](examples/CreateRandomRegions/CreateRandomRegions.md)

    Create a number of regions randomily placed along the genome.
    
     
* [Example 2 - Creating GenomicRanges](examples/CreateGenomicRanges/CreateGenomicRanges.md)

    Create a GenomicRanges (GRanges) object from a data.frame or a local or external file.
    
    
## Complete Examples

* [Example 1 - Multiple data types](Examples/CompleteExamples/MultipleDataTypes/MultipleDataTypes.md)

    A plot with multiple data types plotted together. 
    
    ![Multiple Data Types Figure](Examples/CompleteExamples/MultipleDataTypes/figure/Figure-1.png?raw=true "Multiple Data Types")


*** 

For more extensive documentation on keryoploteR, please refer to the [vignette](http://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf) and the [reference manual](http://bioconductor.org/packages/release/bioc/manuals/karyoploteR/man/karyoploteR.pdf)