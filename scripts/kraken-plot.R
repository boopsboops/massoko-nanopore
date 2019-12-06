#!/usr/bin/env Rscript

library(tidyverse)

sample <- "3mA"
sample <- "22mA"

krep <- read_tsv(file=paste0("../temp-local-only/results/",sample,".out.tsv"),col_names=c("classified","qid","taxid","qlength","lcamap"))

krep %>% ggplot(aes(x=classified,y=qlength,fill=classified)) + geom_boxplot()


ktax <- read_tsv(file=paste0("../temp-local-only/results/",sample,".report.tsv"),col_names=c("percFrag","nClade","nDirect","rank","taxid","scientificName"))


ktax %>% filter(rank=="O") %>% filter(nClade > 100) %>% ggplot(aes(x=fct_reorder(scientificName,nClade,.desc=TRUE),y=nClade,fill=fct_reorder(scientificName,nClade,.desc=TRUE))) + geom_bar(stat="identity")

# combine
ktax %>% mutate(sample=sample) %>% write_csv("../temp-local-only/results/report.combined.tsv",append=TRUE)

# read combined
ktax <- read_csv(file="../temp-local-only/results/report.combined.tsv",col_names=c("percFrag","nClade","nDirect","rank","taxid","scientificName","sample"))

ktax %>% filter(rank=="G") %>% 
    filter(nClade > 200) %>% 
    ggplot(aes(x=fct_reorder(scientificName,nClade,.desc=TRUE),y=nClade,fill=fct_reorder(scientificName,nClade,.desc=TRUE))) + 
    geom_bar(stat="identity") +
    facet_grid(rows=vars(sample))
