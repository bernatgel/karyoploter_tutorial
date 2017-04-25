---
layout: default
title: "karyoploteR Tutorial and Examples"
---

# karyoploteR Tutorial and Examples

This is a tutorial and a few more complex examples on how to use the Bioconductor package karyoploteR to plot linear genomes with arbitrary data.

karryoploteR is an R package to create karyoplots, that is, the chromosome ideograms accompained by any arbitrary data. It is inspired by the R 
base graphics system and does not depend on other graphics packages. 

## Getting Started

karyoploteR is part of [Bioconductor](http://bioconductor.org). The package documentation, including  the [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf)
and [user manual](http://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf) is available at the karyoploteR's 
Bioconductor landing page at [http://bioconductor.org/packages/karyoploteR](http://bioconductor.org/packages/karyoploteR).

To install the package, start R and enter:

{% highlight r %}
  source("https://bioconductor.org/biocLite.R")
  biocLite("karyoploteR")
{% endhighlight %}


## Tutorial

{% for item in site.data.tutorial_and_examples.tutorial %}
    <h3>{{ item.title }}</h3>
      <ul>
        {% for entry in item.subfolderitems %}
          <li><a href="{{ entry.url }}">{{ entry.page }}</a></li>
        {% endfor %}
      </ul>
{% endfor %}


## Complete Examples


