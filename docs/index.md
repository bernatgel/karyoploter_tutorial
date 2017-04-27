---
layout: index
title: "karyoploteR Tutorial and Examples"
---

**[karyoploteR](http://bioconductor.org/packages/karyoploteR)** is an R package to create karyoplots, that is,
representations of whole genomes with arbitrary data plotted on them. It is inspired by the R base graphics system and
does not depend on other graphics packages. The aim of karyoploteR is to offer the user an easy way to plot 
data along the genome to get broad genome-wide view to facilitate the identification of genome wide relations and 
distributions.

**karyoploteR** is based on base R graphics and mimicks its interface. You first create a plot with a call 
to the _plotKaryotype_ function and then sequentially call a number of plotting functions (_kpLines_, _kpPoints_,
_kpBars_…) to add data to the genome plot.

**karyoploteR** is a plotting tool and only a plotting tool. That means that it is not able to download or 
retrieve any data. The downside of this is that the user is responsible of getting the data into R. The upside 
is that it is not tied to any data provider and thus can be used to plot genomic data coming from anywhere.
The only exception to this are the ideograms cytobands, that by default are plotted using predownloaded data
from UCSC.

**karyoploteR** is useful in any situation where a general genome-wide view of data is desirable. It can be
used to plot somatic copy-number changes (SCNA) in cancer genomes obteined from exome, aCGH or SNP-array data;
to plot the global BAM coverage from a WGS experiment; to create manhattan plots from GWAS studies; to create
rainfall plots to detect kataegis. Since it is not tied to any data type or source, karyoploteR can be used to
plot almost anything on a genome-wide scale.


## Getting Started

karyoploteR is part of [Bioconductor](http://bioconductor.org). The package documentation, including  the [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf)
and [user manual](http://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf) is available at the karyoploteR's 
Bioconductor landing page at [http://bioconductor.org/packages/karyoploteR](http://bioconductor.org/packages/karyoploteR).

To install the package, start R and enter:

{% highlight r %}
  source("https://bioconductor.org/biocLite.R")
  biocLite("karyoploteR")
{% endhighlight %}


## <a name="Tutorial"></a>Tutorial

{% for item in site.data.tutorial_and_examples.tutorial %}
  <h3>{{ item.title }}</h3>
  <ul>
    {% for entry in item.subfolderitems %}
      <li><a href="{{ site.baseurl }}/{{ entry.url }}">{{ entry.title }}</a></li>
      <p>{{ entry.text }}</p>
    {% endfor %}
  </ul>
{% endfor %}


## <a name="Examples"></a>Examples


{% for item in site.data.tutorial_and_examples.examples %}
  <div class="col-md-4">
    <div class="thumbnail">
      <a href="{{ site.baseurl }}/{{ item.url }}">
	<img class="img-responsive" src="{{ site.baseurl }}/{{ item.image }}" alt="Image of example {{ item.title }}">
	<div class="caption">
	  <h4>{{ item.title}} </h4>
	  <p>{{ item.text }}</p>
	</div>
      </a>
    </div>
  </div>
{% endfor %}

 
