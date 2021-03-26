# SherryQ
Pathview Code for final visualization of pathways
#https://www.jianshu.com/p/133257d344d3 this is where the intial code is from, but I changed a bit and fixed some bugs to use this for my human dataset. 

#install R packages
# source("https://bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install("gage")
# BiocManager::install("pathview")
# BiocManager::install("gageData")
library("pathview")
library("gage")
library("gageData")
# install.packages("dplyr")
library("dplyr")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(DOSE)
library(stringr)

# BiocManager::install("org.Hsa.eg.db")
library(org.Hs.eg.db)

#data loading
data(kegg.sets.hs)
summary(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs =  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)
head(kegg.sets.hs)
setwd("~/Desktop/RNAseq Pathway analysis 03232021 ")
gene<-read.csv(file="differential_expression_PKD1_vs_control.csv")
head(gene)

gene.df<-bitr(geneID=gene$gene.symbol, fromType = "SYMBOL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(gene.df)
foldchanges = gene$log2FoldChange
names(foldchanges)= gene.df$ENTREZID

head(foldchanges)
keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE,na.rm=TRUE)
keggres=keggres[!is.na(keggres)]
keggres

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

#pick the top and bottom five pathways to draw. you can change the number depends on your interests. 
keggrespathwaysGreater = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tibble::as_tibble() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()

keggrespathwaysLess = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tibble::as_tibble() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()

keggrespathwaysGreater
keggrespathwaysLess

# Get the IDs.
keggresidsGreater = substr(keggrespathwaysGreater, start=1, stop=8)
summary(keggresidsGreater)

keggresidsLess = substr(keggrespathwaysLess, start=1, stop=8)
summary(keggresidsLess)

keggresids=c(keggresidsGreater,keggresidsLess) #prepare to output all of the figures at the sametime

# define the function
plot_pathway = function(pid) 
  pathview(
    gene.data=foldchanges, pathway.id=keggresids, species="hsa", new.signature=FALSE)

# output mutiple pathways and output to the working directory automatically
tmp = sapply(keggresids, plot_pathway ) # takes a bit to finish running, be patient;). 
# note: the human pathway package is a bit outdated, so 
#if it doesnt have the pathways you wanted, you can go to https://pathview.uncc.edu/analysis, somehow this website's library is a bit dated.



#================================================= incase you interested in the other library and cannot install baocManager
# library(go.sets.hs)
# summary(go.sets.hs)
# 
# 
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("StarBioTrek")
# library("StarBioTrek")
# data(kegg.gs)
# data(kegg.gs.dise)
# data(go.gs)
# data(carta.gs)
# 
# species="hsapiens"
# pathwaydb="kegg"
# path<-GetData(species,pathwaydb)

#   install.packages("BiocManager")
