#Iterate through all the directories and knit all Rmd's
#into md's in the docs folder so they can be served by jekyll
#on the GitHub Pages



library(knitr)
library(ezknitr)

base.page.dir <- "../docs"

dirs <- list.dirs(".", recursive = TRUE, full.names = TRUE)

for(d in dirs) {
  if(nchar(d)>2) { #skip . and ..
    message("Working on directory: ", d)
    
    
    files <- list.files(path = d, pattern = "*.Rmd", full.names = TRUE)
    
    for(f in files) {
      message("Knitting file: ", f)
      
      #page.dir <- file.path(base.page.dir, gsub(".Rmd", "", basename(f)))
      page.dir <- file.path(base.page.dir, d)
      dir.create(page.dir, recursive = TRUE)
      
      ezknit(f, out_dir = page.dir, fig_dir="images", keep_html=FALSE)
    }
  }
}

