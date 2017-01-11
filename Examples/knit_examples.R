#Iterate through all the directories and knit all Rmd's


library(knitr)



dirs <- list.dirs(".", recursive = TRUE, full.names = TRUE)

original.wd <- getwd()

for(d in dirs) {
  if(nchar(d)>2) { #skip . and ..
    message("Working on directory: ", d)
    setwd(d)  
    
    files <- list.files(path = ".", pattern = "*.Rmd", full.names = TRUE)
    
    for(f in files) {
      message("Knitting file: ", f)
      knit(f)
    }
    setwd(original.wd)
  }
}
