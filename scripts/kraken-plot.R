#!/usr/bin/env Rscript

# laod packages
library("tidyverse")
library("RColorBrewer")
library("ggthemes")
library("vegan")
library("ggConvexHull")

# load up the results by read
krep <- read_tsv(file="../temp-local-only/results/kraken.results-by-read.tsv",col_names=c("sample","classified","qid","taxid","qlength","lcamap"))

# plot the length of classified and unclassified samples
krep %>% ggplot(aes(x=classified,y=qlength,fill=sample)) + geom_boxplot()

# read the taxonomy summary results
ktax <- read_tsv(file="../temp-local-only/results/kraken.summary.tsv",col_names=c("sample","percFrag","nClade","nDirect","rank","taxid","scientificName"))

# bar plot
ktax %>% 
    mutate(depth=if_else(str_detect(sample,"3"),3,22)) %>% 
    filter(rank=="F") %>% 
    group_by(sample) %>% 
    top_n(n=12,nClade) %>% 
    ggplot(aes(x=fct_reorder(scientificName,nClade,.desc=TRUE),y=nClade,fill=fct_reorder(scientificName,nClade,.desc=TRUE))) + 
    geom_bar(stat="identity") +
    facet_grid(rows=vars(fct_reorder(sample,depth,.desc=FALSE))) #+
    #scale_fill_manual(values=sample(getPalette(n=16)))

# bar with added missing 
ingroup <- ktax %>%
    mutate(rank=str_replace_all(rank,"[0-9]","")) %>%
    filter(rank=="O") %>%
    group_by(sample) %>%
    top_n(n=12,percFrag) %>%
    ungroup() %>%
    select(scientificName) %>%
    distinct() %>%
    pull()

# cols
excols <- colorRampPalette(brewer.pal(11,"Spectral"))(15)

ktax %>%
    group_by(sample) %>%
    mutate(assTot=nClade[rank=="R"],assProp=nClade/assTot) %>%
    ungroup() %>%
    mutate(depth=if_else(str_detect(sample,"3"),3,22)) %>%
    mutate(rank=str_replace_all(rank,"[0-9]","")) %>%
    filter(rank=="O") %>%
    filter(scientificName %in% ingroup) %>%
    group_by(scientificName) %>%
    mutate(taxTot=mean(nClade)) %>%
    ungroup() %>%
    mutate(sample=fct_reorder(sample,depth,.desc=FALSE),scientificName=fct_reorder(scientificName,taxTot,.desc=TRUE)) %>% # reorder the factor
    ggplot(aes(x=scientificName,y=assProp,fill=scientificName)) +
    geom_bar(stat="identity") +
    facet_grid(rows=vars(sample)) +
    scale_fill_manual(values=excols) +
    theme_clean() +
    theme(axis.text.x=element_blank()) +
    labs(fill="Legend",x="Taxonomy",y="Number of reads")


# stacked bar graph
getPalette <- colorRampPalette(brewer.pal(9,"Set1"))
set.seed(1)
ktax %>% 
    mutate(depth=if_else(str_detect(sample,"3"),3,22)) %>% 
    filter(rank=="O") %>%
    group_by(sample) %>%
    top_n(n=30,nClade) %>%
    mutate(nTot=sum(nClade),propAssigned=nClade/nTot) %>%
    ggplot(aes(x=sample,y=propAssigned,fill=fct_reorder(scientificName,propAssigned,.desc=FALSE))) +
    geom_bar(position="fill",stat="identity") + 
    scale_fill_manual(values=sample(getPalette(n=33))) +
    theme_classic()

ktax %>%
    filter(rank=="O") %>%
    #filter(sample!="22mB") %>% 
    #mutate(sampleJoined=if_else(str_detect(sample,"22mB"),"22mB",sample)) %>% 
    #group_by(scientificName,sampleJoined) %>%
    #summarise(nClade=sum(nClade)) %>%
    select(scientificName,sample,nClade) %>% 
    pivot_wider(values_from=nClade,names_from=scientificName,values_fill=list(nClade=0)) %>% 
    data.frame(row.names=1) %>%
    as.matrix() %>% 
    metaMDS(autotransform=TRUE,k=2,distance="bray",binary=FALSE,trace=0) %>%
    vegan::scores() %>%
    as_tibble(rownames="Sample") %>% 
    mutate(depth=as.factor(if_else(str_detect(Sample,"3"),3,22))) %>% 
    ggplot(aes(x=NMDS1,y=NMDS2,color=depth,fill=depth)) + 
    geom_point(size=5,alpha=1,shape=24) + 
    #geom_convexhull(aes(fill=depth),alpha=0.3,colour=NA) +
    geom_text(aes(x=NMDS1,y=NMDS2,color=depth,label=Sample),hjust=1.25,vjust=1.25) +
    scale_fill_ptol() + 
    scale_color_ptol() +
    theme_bw()


# 1. Percentage of fragments covered by the clade rooted at this taxon
# 2. Number of fragments covered by the clade rooted at this taxon
# 3. Number of fragments assigned directly to this taxon
# 4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
#    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
#    Taxa that are not at any of these 10 ranks have a rank code that is
#    formed by using the rank code of the closest ancestor rank with
#    a number indicating the distance from that rank.  E.g., "G2" is a
#    rank code indicating a taxon is between genus and species and the
#    grandparent taxon is at the genus rank.
# 5. NCBI taxonomic ID number
# 6. Indented scientific name
