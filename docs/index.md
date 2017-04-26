---
layout: index
title: "karyoploteR Tutorial and Examples"
---

[karyoploteR](http://bioconductor.org/packages/karyoploteR) is an R package to create karyoplots, that is,
the chromosome ideograms accompained by any arbitrary data. It is inspired by the R base graphics system and
does not depend on other graphics packages. The aim of karyoploteR is to offer the user an easy way to plot 
data along the genome to get broad wide view where it is possible to identify genome wide relations and 
distributions.

karyoploteR is based on base R graphics and mimicks its interface. You first create a plot with plotKaryotype
and then sequentially call a number of functions (kpLines, kpPoints, kpBars…) to add data to the plot.

karyoploteR is a plotting tool and only a plotting tool. That means that it is not able to download or 
retrieve any data. The downside of this is that the user is responsible of getting the data into R. The upside 
is that it is not tied to any data provider and thus can be used to plot genomic data coming from anywhere.
The only exception to this are the ideograms cytobands, that by default are plotted using predownloaded data
from UCSC.

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
      <li><a href="{{ site.baseurl }}/{{ entry.url }}">{{ entry.page }}</a></li>
      <p>{{ entry.text }}</p>
    {% endfor %}
  </ul>
{% endfor %}


## Complete Examples


{% for item in site.data.tutorial_and_examples.examples %}
  <div class="col-md-4">
    <div class="thumbnail">
      <a href="{{ site.baseurl }}/{{ item.url }}">
	<img class="img-responsive" src="{{ item.image }}" alt="Image of example {{ item.title }}">
	<div class="caption">
	  <h4>{{ item.title}} </h4>
	  <p>{{ item.text }}</p>
	</div>
      </a>
    </div>
  </div>
{% endfor %}

 
