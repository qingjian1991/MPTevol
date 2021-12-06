
library(tidyverse)

files = system("ls */*confints_CP.txt", intern = T)
samples = str_split(files, pattern = "/") %>% mapply( function(x) x[1], .)


cp = list()

for(i in 1:length(files)){
  cp[[i]] =  read.table(file = files[i], header = T)[2,] %>%
    mutate(sampleid = samples[i]) 
}

cp = Reduce(rbind, cp)

write.table(cp, file = "summary.cp.txt", quote = F, row.names = F, sep = "\t")

