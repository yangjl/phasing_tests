### Jinliang Yang
### July 30th, 2015


files <- list.files(path="slurm", pattern="txt$")

err <- data.frame()
for(i in 1:length(files)){
    out <- read.table(paste0("test/", files[i]), header=TRUE)
    err <- rbind(err, out)
}

err$er <- err$diff/err$tot


