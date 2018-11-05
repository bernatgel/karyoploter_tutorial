---
layout: index
title: "karyoploteR Tutorial and Examples"
---

**[karyoploteR](http://bioconductor.org/packages/karyoploteR)** is an R package to create karyoplots, that is,
representations of whole genomes with arbitrary data plotted on them. It is inspired by the R base graphics system and
does not depend on other graphics packages. The aim of karyoploteR is to offer the user an easy way to plot 
data along the genome to get broad genome-wide view to facilitate the identification of genome wide relations and 
distributions.

<!-- A carousel showing some of the example images -->
<div id="myCarousel" class="carousel slide" data-interval="false">
 <!-- Indicators -->
 <ol class="carousel-indicators">
  {% assign counter = 1 %}
  {% for item in site.data.tutorial_and_examples.examples %}
   {% if counter == 1 %} 
    <li data-target="#myCarousel" data-slide-to="{{ counter }}" class="active"></li>
   {% else %}
    <li data-target="#myCarousel" data-slide-to="{{ counter }}"></li>
   {% endif %}
   {% assign counter = counter | plus: 1 %}
  {% endfor %}
 </ol>

  <!-- Wrapper for slides -->
 <div class="carousel-inner">
  {% assign counter = 1 %}
  {% for item in site.data.tutorial_and_examples.examples %}
    {% if counter == 1 %}
      <div class="item active">
    {% else %}
      <div class="item">
    {% endif %}
    {% assign counter = counter | plus: 1 %}
     <img class="carousel-img" src="{{ site.baseurl }}/{{ item.image }}" alt="Image of example {{ item.title }}">
     {% if item.only_devel == 1 %}
	<div class="devel-only">devel only</div>
     {% endif %}
     <div class="carousel-caption">
      <h4>Example - {{ item.title}}</h4>
      <a class="btn btn-lg btn-primary" href="{{ site.baseurl }}/{{ item.url }}" role="button">See Example</a>
     </div>
    </div>
  {% endfor %}
  </div>

  <!-- Left and right controls -->
  <a class="left carousel-control" href="#myCarousel" data-slide="prev">
    <span class="glyphicon glyphicon-chevron-left"></span>
    <span class="sr-only">Previous</span>
  </a>
  <a class="right carousel-control" href="#myCarousel" data-slide="next">
    <span class="glyphicon glyphicon-chevron-right"></span>
    <span class="sr-only">Next</span>
  </a>
</div>

&nbsp;

**karyoploteR** is based on base R graphics and mimicks its interface. You first create a plot with a call 
to the `plotKaryotype` function and then sequentially call a number of plotting functions (`kpLines`, `kpPoints`,
`kpBars`…) to add data to the genome plot.

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


## <a name="GettingStarted"></a>Getting Started

karyoploteR is part of [Bioconductor](http://bioconductor.org) since version BioC 3.5. The package documentation, including  the [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf)
and [user manual](http://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf) is available at the karyoploteR's 
Bioconductor landing page at [http://bioconductor.org/packages/karyoploteR](http://bioconductor.org/packages/karyoploteR).

To install the package you'll need to use [Bioconductor's own package manager](https://www.bioconductor.org/install/), called [`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/index.html).

To do so, simply start R and enter the following code:

{% highlight r %}
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("karyoploteR")
{% endhighlight %}



#### Usign the development version

To use the [development version of karyoploteR](http://bioconductor.org/packages/devel/bioc/html/karyoploteR.html) 
you should use the [devel version of Bioconductor](https://www.bioconductor.org/developers/how-to/useDevel/). The 
devel version of the package might work with release version of Bioconductor, althought that's not expected to be
always the case. You should be able to install the development version from the 
[github repo](https://github.com/bernatgel/karyoploter) using `install_github()`
from the [devtools](https://github.com/hadley/devtools) package.

## <a name="Citing"></a>Citing karyoploteR

karyoploteR has been developed by **Bernat Gel** [<i class="fa fa-envelope" aria-hidden="true"></i>](mailto:bgel@igtp.cat)
[<i class="fa fa-twitter" aria-hidden="true"></i>](https://twitter.com/bernatgel) 
[<i class="fa fa-github" aria-hidden="true"></i>](https://github.com/bernatgel)
  and   **Eduard Serra** [<i class="fa fa-envelope" aria-hidden="true"></i>](mailto:eserra@igtp.cat) at [IGTP](http://www.germanstrias.org/)
Hereditary Cancer Group.

If you use karyoploteR in your research, please cite the Bioinformatics paper describing it:

Bernat Gel & Eduard Serra. (2017). *karyoploteR: an R/Bioconductor package to plot customizable genomes displaying arbitrary data*. Bioinformatics, 31–33. [doi:10.1093/bioinformatics/btx346](https://doi.org/10.1093/bioinformatics/btx346)




## <a name="Tutorial"></a>Tutorial

The tutorial is a work in progress yet. Feel free to [contact us](mailto:bgel@igtp.cat) to ask for any clarification or propose a a new section.

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
      {% if item.only_devel == 1 %}
	<div class="devel-only">devel only</div>
      {% endif %}
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



 
