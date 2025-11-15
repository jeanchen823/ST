# SpaGene

```bash
https://genome.cshlp.org/content/32/9/1736

https://github.com/liuqivandy/SpaGene

library(devtools)
install_github("liuqivandy/SpaGene")

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(SeuratObject)
library(SpaGene)
library(tidyverse) 

setwd("BreastCancer/")
load("bc_raw.rds")

bc_spagene<-SpaGene(count,location)
# the most signifiant spatially variable genes
head(bc_spagene$spagene_res[order(bc_spagene$spagene_res$adjp),])

pattern<-FindPattern(bc_spagene)
PlotPattern(pattern,location)

top5<-apply(pattern$genepattern,2,function(x){names(x)[order(x,decreasing=T)][1:5]})
library(pheatmap)
pheatmap(pattern$genepattern[rownames(pattern$genepattern)%in%top5,])

load("LRpair_human.rds")
bc_lr<-SpaGene_LR(count,location,LRpair=LRpair) 
# the most signficant colocalized LR pairs
head(bc_lr[order(bc_lr$adj),])


plotLR(count,location,LRpair=c("FN1","SDC2"),alpha.min=0.5)

https://www.nature.com/articles/s41467-024-51580-7


https://github.com/JiaLiVUMC/GCA_LND/blob/main/Manuscript/spatial.R

```



```bash
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(SpaGene)
#Seurat 
library(Seurat)
library(SeuratObject)
library(SpaGene)
library(tidyverse) 
setwd("brain10X/")
load("brain10x_raw.rds")

brain10x_sv<-SpaGene(count,location)

pattern<-FindPattern(brain10x_sv,nPattern = 15)

PlotPattern(pattern,location,pt.size = 0.5)

top5<-apply(pattern$genepattern,2,function(x){names(x)[order(x,decreasing=T)][1:5]})
library(pheatmap)
pheatmap(pattern$genepattern[rownames(pattern$genepattern)%in%top5,],fontsize_row = 6)


library(SpaGene)
library(SeuratData)
library(Seurat)


brain1 <- LoadData("stxBrain", type = "anterior1")
save(brain1,file = "brain1.rdata")
brain2 <- LoadData("stxBrain", type = "posterior1")
save(brain2,file = "brain2.rdata")
count1<-GetAssayData(brain1,slot="counts")
count2<-GetAssayData(brain2,slot="counts")


location1<-GetTissueCoordinates(brain1)
location2<-GetTissueCoordinates(brain2)

spa1<-SpaGene(count1,location1)
spa2<-SpaGene(count2,location2)

pattern<-FindPattern_Multi(list(spa1,spa2),nPattern=25)

locationlist<-list(location1[,2:1],location2[,2:1])
patternnum<-nrow(pattern$pattern)
for (i in 1:10) {
  cat(paste0("pattern ",i,"\n"))
  
  print(PlotPattern_Multi(pattern,locationlist,pt.size=1,patternid=i,max.cutoff =0.95))
}




```



```bash
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())


library(Seurat)
library(SeuratObject)
library(SpaGene)
library(tidyverse) 

setwd("BreastCancer/")
load("bc_raw.rds")

bc_spagene<-SpaGene(count,location)
# the most signifiant spatially variable genes
head(bc_spagene$spagene_res[order(bc_spagene$spagene_res$adjp),])

pattern<-FindPattern(bc_spagene)

PlotPattern(pattern,location)

top5<-apply(pattern$genepattern,2,function(x){names(x)[order(x,decreasing=T)][1:5]})
library(pheatmap)
pheatmap(pattern$genepattern[rownames(pattern$genepattern)%in%top5,])


load("LRpair_human.rds")
bc_lr<-SpaGene_LR(count,location,LRpair=LRpair) 
# the most signficant colocalized LR pairs
head(bc_lr[order(bc_lr$adj),])


plotLR(count,location,LRpair=c("FN1","SDC2"),alpha.min=0.5)
```



```bash
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())


library(Seurat)
library(SeuratObject)
library(SpaGene)
library(tidyverse) 


load("brain10x_raw.rds")

### Find spatially variable genes and patterns
brain10x_sv<-SpaGene(count,location)
# the most significant spt
head(brain10x_sv$spagene_res[order(brain10x_sv$spagene_res$adjp),])

pattern<-FindPattern(brain10x_sv,nPattern = 15)

PlotPattern(pattern,location,pt.size = 0.5)

## Top 5 genes falling into each pattern
top5<-apply(pattern$genepattern,2,function(x){names(x)[order(x,decreasing=T)][1:5]})
library(pheatmap)
pheatmap(pattern$genepattern[rownames(pattern$genepattern)%in%top5,],fontsize_row = 6)



library(SpaGene)
library(SeuratData)
library(Seurat)

brain1 <- LoadData("stxBrain", type = "anterior1")
save(brain1,file = "brain1.rdata")
brain2 <- LoadData("stxBrain", type = "posterior1")
save(brain2,file = "brain2.rdata")
count1<-GetAssayData(brain1,slot="counts")
count2<-GetAssayData(brain2,slot="counts")


location1<-GetTissueCoordinates(brain1)
location2<-GetTissueCoordinates(brain2)

spa1<-SpaGene(count1,location1)
spa2<-SpaGene(count2,location2)

pattern<-FindPattern_Multi(list(spa1,spa2),nPattern=25)


locationlist<-list(location1[,2:1],location2[,2:1])
patternnum<-nrow(pattern$pattern)
for (i in 1:10) {
  cat(paste0("pattern ",i,"\n"))
  
  print(PlotPattern_Multi(pattern,locationlist,pt.size=1,patternid=i,max.cutoff =0.95))
}



```



