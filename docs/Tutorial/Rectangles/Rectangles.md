---
layout: tutorial
label: Rectangles
title: Rectangles
---




## Plotting Rectangles

To plot rectangles in a karyoplot we need to use the `kpRect` function. Given a
character vector _chr_ and 4 integer vectors _x0_, _y0_, _x1_ and _y1_ it 
will plot rectangles with an angle at (_x0_, _y0_) and the oposite one at 
(_x1_, _y1_).


```r
library(karyoploteR)

kp <- plotKaryotype(chromosomes="chr1")
kpRect(kp, chr="chr1", x0=0, x1=80e6, y0=0.2, y1=0.8)
```

![plot of chunk Figure1](images//Figure1-1.png)

We can give it vectors of positions and it will plot a rectangle for each element 
in the vectors (recycling them if necessary).


```r
x0 <- 1:23*10e6
x1 <- 2:24*10e6
y0 <- rnorm(23, mean=0.3, sd=0.1)
y1 <- c(0.7, 0.9)

kp <- plotKaryotype(chromosomes="chr1")
kpRect(kp, chr="chr1", x0=x0, x1=x1, y0=y0, y1=y1)
```

![plot of chunk Figure2](images//Figure2-1.png)

The rectangles can be customized with the same 
[graphical parameters](https://www.rdocumentation.org/packages/graphics/topics/par)
as in the R base graphics `rect` function: _border_, _lty_, _col_... Taking into
account that to produce rectangles without border we need `border=NA` (and not
`border=NULL`).


```r
kp <- plotKaryotype(chromosomes="chr1")
kpRect(kp, chr="chr1", x0=x0, x1=x1, y0=y0, y1=y1, col=rainbow(23), 
           lty=c(1,2,3,4), lwd=(1:23)/4, r0=0.5, r1=1, border=c("black", NA))
kpRect(kp, chr="chr1", x0=x0, x1=x1, y0=y0, y1=y1, border=rainbow(23), 
           lty=c(1,2,3,4), lwd=(1:23)/4, r0=0.5, r1=0)
```

![plot of chunk Figure3](images//Figure3-1.png)





