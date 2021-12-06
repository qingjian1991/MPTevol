library(sequenza)

library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("--seqz.file"),  default=NULL ,type = "character", help = "seqz file. [Required]"),
  make_option(c("--sample.id"),  default=NULL ,type = "character", help = "sample id. [Required]")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

#opt = list(seqz.file = "Breast_1.bin50.seqz.gz", sample.id = "Breast_1" )

seqz.file = opt$seqz.file
sample.id = opt$sample.id

seqz <- sequenza.extract(seqz.file, verbose = FALSE, parallel = 20)


seqz.fit <- sequenza.fit(seqz, mc.cores = 20)

sequenza.results(sequenza.extract = seqz,
                 cp.table = seqz.fit, sample.id = sample.id,
                 out.dir= sample.id)


