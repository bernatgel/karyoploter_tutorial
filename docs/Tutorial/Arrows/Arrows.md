---
layout: tutorial
label: Arrows
title: Arrows
---




## Plotting Arrows

To plot arrows in a karyoplot we need to use the `kpArrows` function. Given a
character vector _chr_ and 4 integer vectors _x0_, _y0_, _x1_ and _y1_ it 
will plot arrows going from (_x0_, _y0_) to (_x1_, _y1_).


```r
library(karyoploteR)

kp <- plotKaryotype(chromosomes="chr1")
kpArrows(kp, chr="chr1", x0=0, x1=80e6, y0=0.2, y1=0.8)
```

![plot of chunk Figure1](images//Figure1-1.png)

We can give it vectors of positions and it will plot an arrow for each element 
in the vectors (recycling them if necessary).


```r
x0 <- 1:23*10e6
x1 <- 2:24*10e6
y0 <- rnorm(23, mean=0.3, sd=0.1)
y1 <- c(0.7, 0.9)

kp <- plotKaryotype(chromosomes="chr1")
kpArrows(kp, chr="chr1", x0=x0, x1=x1, y0=y0, y1=y1)
```

![plot of chunk Figure2](images//Figure2-1.png)

The lines can be customized with the same 
[graphical parameters](https://www.rdocumentation.org/packages/graphics/topics/par)
as in the R base graphics `arrows` function: _lwd_, _lty_, _col_, _length_,
_angle_...


```r
kp <- plotKaryotype(chromosomes="chr1")
kpArrows(kp, chr="chr1", x0=x0, x1=x1, y0=0.5, y1=y1, col=rainbow(23), 
          lwd=(1:23)/6, angle=30, length=0.2)
kpArrows(kp, chr="chr1", x0=x0, x1=x0, y0=0.5, y1=0, col=rainbow(23), 
          lwd=(1:23)/6, angle=10, length=0.1)
```

![plot of chunk Figure3](images//Figure3-1.png)





