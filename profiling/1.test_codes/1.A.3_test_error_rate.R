### Jinliang Yang
### July 30th, 2015


files <- list.files(path="largedata/out/", pattern="txt$")

out <- data.frame()
for(i in 1:length(files)){
    tem <- read.table(paste0("largedata/out/", files[i]), header=TRUE)
    tem$file <- files[i]
    out <- rbind(out, tem)
}

err$er <- err$err/err$tot


