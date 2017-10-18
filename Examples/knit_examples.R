#Iterate through all the directories and knit all Rmd's
#into md's in the docs folder so they can be served by jekyll
#on the GitHub Pages

library(pasilla)

library(knitr)
library(ezknitr)

force.knit.all <- FALSE



base.page.dir <- "../docs"

dirs <- list.dirs(".", recursive = TRUE, full.names = TRUE)

for(d in dirs) {
  if(nchar(d)>2) { #skip . and ..
    message("Working on directory: ", d)
    
    
    files <- list.files(path = d, pattern = "*.Rmd", full.names = TRUE)
    
    for(f in files) {
      message("Checking source file: ", f)
      
      #page.dir <- file.path(base.page.dir, gsub(".Rmd", "", basename(f)))
      page.dir <- file.path(base.page.dir, d)
      
      if(dir.exists(page.dir)) {
        dest.file <- file.path(page.dir, gsub("Rmd", "md", basename(f)))
        if(file.exists(dest.file)) {
          #If the creation of the knitted file is posterior to the modification
          #of the source file, do not knit
          if(file.info(dest.file)$ctime > file.info(f)$mtime & !force.knit.all) {
            message("No changes in ", f, ". Skipping.")
            next;
          }
       } 
        
      } else {
        dir.create(page.dir, recursive = TRUE)
      }
      
      message("Knitting ", f)
      ezknit(f, out_dir = page.dir, fig_dir="images", keep_html=FALSE)
    }
  }
}

