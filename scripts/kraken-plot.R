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


# choose rank
srank <- "O"
topn <- 12

# get all data from top n, including missing 
ingroup <- ktax %>%
    mutate(rank=str_replace_all(rank,"[0-9]","")) %>%
    filter(rank==srank) %>%
    group_by(sample) %>%
    top_n(n=topn,percFrag) %>%
    ungroup() %>%
    select(scientificName) %>%
    distinct() %>%
    pull()

# get some new colours
excols <- colorRampPalette(brewer.pal(11,"Spectral"))(17)

# run the barplot of proportions
ktax %>%
    group_by(sample) %>%
    mutate(assTot=nClade[rank=="R"],assProp=nClade/assTot) %>%
    ungroup() %>%
    mutate(depth=if_else(str_detect(sample,"3"),3,22)) %>%
    mutate(rank=str_replace_all(rank,"[0-9]","")) %>%
    filter(rank==srank) %>%
    filter(scientificName %in% ingroup) %>%
    group_by(scientificName) %>%
    mutate(taxTot=mean(nClade)) %>%
    ungroup() %>%
    mutate(sample=fct_reorder(sample,depth,.desc=FALSE),scientificName=fct_reorder(scientificName,taxTot,.desc=TRUE)) %>%
    ggplot(aes(x=scientificName,y=assProp,fill=scientificName)) +
    geom_bar(stat="identity") +
    facet_grid(rows=vars(sample)) +
    scale_fill_manual(values=excols) +
    theme_light() +
    theme(axis.text.x=element_blank(),panel.grid.major.x=element_blank()) +
    labs(fill="Legend",x="Taxonomy",y="Proportion of assigned reads")


# make a stacked bar graph - not great
getPalette <- colorRampPalette(brewer.pal(9,"Set1"))
set.seed(1)
ktax %>% 
    mutate(depth=if_else(str_detect(sample,"3"),3,22)) %>% 
    filter(rank==srank) %>%
    group_by(sample) %>%
    top_n(n=topn,nClade) %>%
    mutate(nTot=sum(nClade),propAssigned=nClade/nTot) %>%
    ggplot(aes(x=sample,y=propAssigned,fill=fct_reorder(scientificName,propAssigned,.desc=FALSE))) +
    geom_bar(position="fill",stat="identity") + 
    scale_fill_manual(values=sample(getPalette(n=33))) +
    theme_classic()

# NMDS
ktax %>%
    group_by(sample) %>%
    mutate(assTot=nClade[rank=="R"],assProp=nClade/assTot) %>%
    ungroup() %>%
    mutate(rank=str_replace_all(rank,"[0-9]","")) %>%
    filter(rank==srank) %>%
    #filter(sample!="22mB") %>% 
    #mutate(sampleJoined=if_else(str_detect(sample,"22mB"),"22mB",sample)) %>% 
    #group_by(scientificName,sampleJoined) %>%
    #summarise(nClade=sum(nClade)) %>%
    select(scientificName,sample,assProp) %>% 
    pivot_wider(values_from=assProp,names_from=scientificName,values_fill=list(assProp=0)) %>% 
    data.frame(row.names=1) %>%
    as.matrix() %>% 
    metaMDS(autotransform=FALSE,k=2,distance="bray",binary=FALSE,trace=2) %>%
    vegan::scores() %>%
    as_tibble(rownames="Sample") %>% 
    mutate(depth=as.factor(if_else(str_detect(Sample,"3"),3,22))) %>% 
    ggplot(aes(x=NMDS1,y=NMDS2,color=depth,fill=depth)) + 
    geom_point(size=5,alpha=1,shape=24) + 
    #geom_convexhull(aes(fill=depth),alpha=0.3,colour=NA) +
    geom_text(aes(x=NMDS1,y=NMDS2,color=depth,label=Sample),position=position_nudge(y=-0.004),check_overlap=FALSE) +#hjust=-0.025,vjust=1.75
    scale_fill_ptol() + 
    scale_color_ptol() +
    theme_bw()


# summary stats
ktax.tab <- ktax %>%
    #mutate(rank=str_replace_all(rank,"[0-9]","")) %>%
    filter(str_detect(rank,"[0-9]",negate=TRUE)) %>% 
    group_by(sample,rank) %>%
    summarise(n=sum(nClade)) %>%
    ungroup() %>% 
    pivot_wider(values_from=n,names_from=rank) %>% 
    mutate(T=U+R) %>% # just to check
    pivot_longer(cols=-sample,names_to="rank",values_to="n") %>% 
    mutate(rank=fct_relevel(rank,c("T","U","R","D","K","P","C","O","F","G","S"))) %>%
    arrange(sample,rank) %>%
    pivot_wider(values_from=n,names_from=sample) %>% 
    filter(rank!="K") %>%
    select(rank,`3mA`,`3mB`,`22mA`,`22mB`,`22mBv2`) %>%
    filter(rank!="T" & rank!="U" & rank!="R") # now ditch

# summary
krep.tab <- krep %>%
    group_by(sample) %>%
    summarise(nFilteredReads=length(classified),maxLengthFiltered=max(qlength),meanLengthFiltered=mean(qlength),medianLengthFiltered=median(qlength),unclassified=length(which(classified=="U")),classified=length(which(classified=="C"))) %>%
    mutate(nRawReads=c(155906,118847,186798,199847,170166),filteredN50=c(3122,1526,2069,5289,8566)) %>%
    select(sample,nRawReads,nFilteredReads,filteredN50,meanLengthFiltered,medianLengthFiltered,maxLengthFiltered,unclassified,classified) %>%
    pivot_longer(cols=-sample,names_to="rank",values_to="n") %>% 
    pivot_wider(values_from=n,names_from=sample) %>%
    select(rank,`3mA`,`3mB`,`22mA`,`22mB`,`22mBv2`)

# join the two tables (ignore warning)
tab.combined <- bind_rows(krep.tab,ktax.tab)

# write out
#write_csv(tab,path="../temp-local-only/summary-table.csv")
write_csv(tab.combined,path="figures/summary-table.csv")


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
