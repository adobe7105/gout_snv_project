
library(dplyr)
library(readxl)
setwd("/run/media/r730/KESU/instrain_data")
#1测试组test基因上snv数量######################################################

path      <- "/run/media/r730/KESU/instrain_data/test_output_gene_control"
fileNames <- dir(path)
filePath  <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})  
test_output_gene_control <- lapply(filePath, function(x){
  read.delim(x)})

path2      <- "/run/media/r730/KESU/instrain_data/test_output_gene_gout"
fileNames2 <- dir(path2)
filePath2  <- sapply(fileNames2, function(x){ 
  paste(path2,x,sep='/')})  
test_output_gene_gout <- lapply(filePath2, function(x){
  read.delim(x)})
#2训练组find基因上snv数量######################################################

path3      <- "/run/media/r730/KESU/instrain_data/find_output_gene_0w_control"
fileNames3 <- dir(path3)
filePath3  <- sapply(fileNames3, function(x){ 
  paste(path3,x,sep='/')})  
find_output_gene_0w_control <- lapply(filePath3, function(x){
  read.delim(x)})

path4      <- "/run/media/r730/KESU/instrain_data/find_output_gene_0w_gout"
fileNames4 <- dir(path4)
filePath4  <- sapply(fileNames4, function(x){ 
  paste(path4,x,sep='/')})  
find_output_gene_0w_gout <- lapply(filePath4, function(x){
  read.delim(x)})
#3测试组test基因组X10######################################################

path5      <- "/run/media/r730/KESU/instrain_data/test_output_genome_control"
fileNames5 <- dir(path5)
filePath5  <- sapply(fileNames5, function(x){ 
  paste(path5,x,sep='/')})  
test_output_genome_control <- lapply(filePath5, function(x){
  read.delim(x)})

path6      <- "/run/media/r730/KESU/instrain_data/test_output_genome_gout"
fileNames6 <- dir(path6)
filePath6  <- sapply(fileNames6, function(x){ 
  paste(path6,x,sep='/')})  
test_output_genome_gout <- lapply(filePath6, function(x){
  read.delim(x)})
#4训练组find基因组X10######################################################

path7      <- "/run/media/r730/KESU/instrain_data/find_output_genome_0w_control"
fileNames7 <- dir(path7)
filePath7  <- sapply(fileNames7, function(x){ 
  paste(path7,x,sep='/')})  
find_output_genome_control <- lapply(filePath7, function(x){
  read.delim(x)})

path8      <- "/run/media/r730/KESU/instrain_data/find_output_genome_0w_gout"
fileNames8 <- dir(path8)
filePath8  <- sapply(fileNames8, function(x){ 
  paste(path8,x,sep='/')})  
find_output_genome_gout <- lapply(filePath8, function(x){
  read.delim(x)})

#5训练组snv位置######################################################

path9      <- "/run/media/r730/KESU/instrain_data/find_output_snv_0w_control"
fileNames9 <- dir(path9)
filePath9  <- sapply(fileNames9, function(x){ 
  paste(path9,x,sep='/')})  
find_output_snv_control <- lapply(filePath9, function(x){
  read.delim(x)})

path10      <- "/run/media/r730/KESU/instrain_data/find_output_snv_0w_gout"
fileNames10 <- dir(path10)
filePath10  <- sapply(fileNames10, function(x){ 
  paste(path10,x,sep='/')})  
find_output_snv_gout <- lapply(filePath10, function(x){
  read.delim(x)})
#6测试组snv位置######################################################

path11      <- "/run/media/r730/KESU/instrain_data/test_output_snv_control"
fileNames11 <- dir(path11)
filePath11  <- sapply(fileNames11, function(x){ 
  paste(path11,x,sep='/')})  
test_output_snv_control <- lapply(filePath11, function(x){
  read.delim(x)})

path12      <- "/run/media/r730/KESU/instrain_data/test_output_snv_gout"
fileNames12 <- dir(path12)
filePath12  <- sapply(fileNames12, function(x){ 
  paste(path12,x,sep='/')})  
test_output_snv_gout <- lapply(filePath12, function(x){
  read.delim(x)})
#保存rds格式###############################################
saveRDS(find_output_genome_gout,"find_output_genome_gout.rds")
saveRDS(find_output_genome_control,"find_output_genome_control.rds")
saveRDS(test_output_genome_gout,"test_output_genome_gout.rds")
saveRDS(test_output_genome_control,"test_output_genome_control.rds")

saveRDS(find_output_gene_0w_gout,"find_output_gene_0w_gout.rds")
saveRDS(find_output_gene_0w_control,"find_output_gene_0w_control.rds")
saveRDS(test_output_gene_gout,"test_output_gene_gout.rds")
saveRDS(test_output_gene_control,"test_output_gene_control.rds")

saveRDS(find_output_snv_gout,"find_output_snv_gout.rds")
saveRDS(find_output_snv_control,"find_output_snv_control.rds")
saveRDS(test_output_snv_gout,"test_output_snv_gout.rds")
saveRDS(test_output_snv_control,"test_output_snv_control.rds")
#读取rds格式##############################################
library(dplyr)
library(readxl)
setwd("/run/media/r730/KESU/instrain_data")

find_output_genome_gout<-readRDS("find_output_genome_gout.rds")
find_output_genome_control<-readRDS("find_output_genome_control.rds")
test_output_genome_gout<-readRDS("test_output_genome_gout.rds")
test_output_genome_control<-readRDS("test_output_genome_control.rds")

find_output_gene_0w_gout<-readRDS("find_output_gene_0w_gout.rds")
find_output_gene_0w_control<-readRDS("find_output_gene_0w_control.rds")
test_output_gene_gout<-readRDS("test_output_gene_gout.rds")
test_output_gene_control<-readRDS("test_output_gene_control.rds")


find_output_snv_gout<- readRDS("find_output_snv_gout.rds")
find_output_snv_control<- readRDS("find_output_snv_control.rds")
test_output_snv_gout<- readRDS("test_output_snv_gout.rds")
test_output_snv_control<-readRDS("test_output_snv_control.rds")

#gene合并对应的coverage################################################

genomes <- read.delim("genomes.stb", header=FALSE)
colnames(genomes)=c("scaffold","genome")
species <- read.delim("species.txt")
nnn <- read_excel("nnn.xlsx",sheet = "genomes")
colnames(nnn)=c("scaffold","clade_name")
genome_scaffold_name<-left_join(genomes,nnn,by="scaffold")

name_77=c()
name_63=c()
name_23=c()
name_25=c()

#整理合并，选择是否进行调整####################################################################

for (i in 1:77) {
  
  name_s=strsplit(names(find_output_gene_0w_gout[i]),split = ".Re")[[1]][1]
  fengdu_species_clade=species[,c("clade_name",name_s)]
  genome_name=find_output_genome_gout[[i]][,c("genome","coverage")]
  genome_fengdu=species[,c("clade_name",name_s)]
  colnames(genome_fengdu)=c("clade_name","dd")
  genome_coverage_fengdu1=data.frame()
  genome_coverage_fengdu2=data.frame()
  genome_coverage_fengdu1=left_join(genome_scaffold_name,genome_fengdu,by="clade_name")
  genome_coverage_fengdu2=left_join(genome_coverage_fengdu1,genome_name,by="genome")
  find_output_gene_0w_gout[[i]]=left_join(find_output_gene_0w_gout[[i]],genome_coverage_fengdu2,by="scaffold")
  find_output_gene_0w_gout[[i]]=find_output_gene_0w_gout[[i]]%>%filter(coverage.y>=10)
  find_output_gene_0w_gout[[i]]=find_output_gene_0w_gout[[i]]%>%filter(dd>0)
  find_output_gene_0w_gout[[i]]$snv_all=(find_output_gene_0w_gout[[i]]$SNV_count+find_output_gene_0w_gout[[i]]$SNS_count)
  find_output_gene_0w_gout[[i]]$standard_snv_all=(find_output_gene_0w_gout[[i]]$SNV_count+find_output_gene_0w_gout[[i]]$SNS_count)/
    (find_output_gene_0w_gout[[i]]$dd*find_output_gene_0w_gout[[i]]$coverage.y)
  find_output_gene_0w_gout[[i]]=find_output_gene_0w_gout[[i]][,c("gene","standard_snv_all")]
  name_77=c(name_77,name_s)
}


for (i in 1:63) {
  
  name_s=strsplit(names(find_output_gene_0w_control[i]),split = ".Re")[[1]][1]
  fengdu_species_clade=species[,c("clade_name",name_s)]
  genome_name=find_output_genome_control[[i]][,c("genome","coverage")]
  genome_fengdu=species[,c("clade_name",name_s)]
  colnames(genome_fengdu)=c("clade_name","dd")
  genome_coverage_fengdu1=data.frame()
  genome_coverage_fengdu2=data.frame()
  genome_coverage_fengdu1=left_join(genome_scaffold_name,genome_fengdu,by="clade_name")
  genome_coverage_fengdu2=left_join(genome_coverage_fengdu1,genome_name,by="genome")
  find_output_gene_0w_control[[i]]=left_join(find_output_gene_0w_control[[i]],genome_coverage_fengdu2,by="scaffold")
  find_output_gene_0w_control[[i]]=find_output_gene_0w_control[[i]]%>%filter(coverage.y>=10)
  find_output_gene_0w_control[[i]]=find_output_gene_0w_control[[i]]%>%filter(dd>0)
  find_output_gene_0w_control[[i]]$snv_all=(find_output_gene_0w_control[[i]]$SNV_count+find_output_gene_0w_control[[i]]$SNS_count)
  find_output_gene_0w_control[[i]]$standard_snv_all=(find_output_gene_0w_control[[i]]$SNV_count+find_output_gene_0w_control[[i]]$SNS_count)/
    (find_output_gene_0w_control[[i]]$dd*find_output_gene_0w_control[[i]]$coverage.y)
  find_output_gene_0w_control[[i]]=find_output_gene_0w_control[[i]][,c("gene","standard_snv_all")]
  name_63=c(name_63,name_s)
}

for (i in 1:23) {
  
  name_s=strsplit(names(test_output_gene_control[i]),split = ".Re")[[1]][1]
  fengdu_species_clade=species[,c("clade_name",name_s)]
  genome_name=test_output_genome_control[[i]][,c("genome","coverage")]
  genome_fengdu=species[,c("clade_name",name_s)]
  colnames(genome_fengdu)=c("clade_name","dd")
  genome_coverage_fengdu1=data.frame()
  genome_coverage_fengdu2=data.frame()
  genome_coverage_fengdu1=left_join(genome_scaffold_name,genome_fengdu,by="clade_name")
  genome_coverage_fengdu2=left_join(genome_coverage_fengdu1,genome_name,by="genome")
  test_output_gene_control[[i]]=left_join(test_output_gene_control[[i]],genome_coverage_fengdu2,by="scaffold")
  test_output_gene_control[[i]]=test_output_gene_control[[i]]%>%filter(coverage.y>=10)
  test_output_gene_control[[i]]=test_output_gene_control[[i]]%>%filter(dd>0)
  test_output_gene_control[[i]]$snv_all=(test_output_gene_control[[i]]$SNV_count+test_output_gene_control[[i]]$SNS_count)
  test_output_gene_control[[i]]$standard_snv_all=(test_output_gene_control[[i]]$SNV_count+test_output_gene_control[[i]]$SNS_count)/
    (test_output_gene_control[[i]]$dd*test_output_gene_control[[i]]$coverage.y)
  test_output_gene_control[[i]]=test_output_gene_control[[i]][,c("gene","standard_snv_all")]
  name_23=c(name_23,name_s)
}


for (i in 1:25) {
  
  name_s=strsplit(names(test_output_gene_gout[i]),split = ".Re")[[1]][1]
  fengdu_species_clade=species[,c("clade_name",name_s)]
  genome_name=test_output_genome_gout[[i]][,c("genome","coverage")]
  genome_fengdu=species[,c("clade_name",name_s)]
  colnames(genome_fengdu)=c("clade_name","dd")
  genome_coverage_fengdu1=data.frame()
  genome_coverage_fengdu2=data.frame()
  genome_coverage_fengdu1=left_join(genome_scaffold_name,genome_fengdu,by="clade_name")
  genome_coverage_fengdu2=left_join(genome_coverage_fengdu1,genome_name,by="genome")
  test_output_gene_gout[[i]]=left_join(test_output_gene_gout[[i]],genome_coverage_fengdu2,by="scaffold")
  test_output_gene_gout[[i]]=test_output_gene_gout[[i]]%>%filter(coverage.y>=10)
  test_output_gene_gout[[i]]=test_output_gene_gout[[i]]%>%filter(dd>0)
  test_output_gene_gout[[i]]$snv_all=(test_output_gene_gout[[i]]$SNV_count+test_output_gene_gout[[i]]$SNS_count)
  test_output_gene_gout[[i]]$standard_snv_all=(test_output_gene_gout[[i]]$SNV_count+test_output_gene_gout[[i]]$SNS_count)/
    (test_output_gene_gout[[i]]$dd*test_output_gene_gout[[i]]$coverage.y)
  test_output_gene_gout[[i]]=test_output_gene_gout[[i]][,c("gene","standard_snv_all")]
  name_25=c(name_25,name_s)
}

for (i in 2:77) { 
  find_output_gene_0w_gout[[1]] <-full_join(find_output_gene_0w_gout[[1]], find_output_gene_0w_gout[[i]], by = "gene")
}

for (i in 2:63) { 
  find_output_gene_0w_control[[1]] <- full_join(find_output_gene_0w_control[[1]], find_output_gene_0w_control[[i]], by = "gene")
}

for (i in 2:25) { 
  test_output_gene_gout[[1]] <- full_join(test_output_gene_gout[[1]], test_output_gene_gout[[i]], by = "gene")
}

for (i in 2:23) { 
  test_output_gene_control[[1]] <- full_join(test_output_gene_control[[1]], test_output_gene_control[[i]], by = "gene")
}

data77=find_output_gene_0w_gout[[1]]
data63=find_output_gene_0w_control[[1]]
data23=test_output_gene_control[[1]]
data25=test_output_gene_gout[[1]]

for (i in 2:64) { 
  names(data63)[i]<-name_63[i-1]}
for (i in 2:78) { 
  names(data77)[i]<-name_77[i-1]}
for (i in 2:24) { 
  names(data23)[i]<-name_23[i-1]}
for (i in 2:26) { 
  names(data25)[i]<-name_25[i-1]}

data_all48=full_join(data23,data25,by="gene")
data_all140=full_join(data63,data77,by="gene")
data_188=full_join(data_all48,data_all140,by="gene")

#处理数量数据框#######################################################################

data_find=full_join(data63,data77,by="gene")
data_find[is.na(data_find)]<-0
#data_find=na.omit(data_find)
rownames(data_find)=data_find[,1]
data_find=data_find[,-1]
data_find=t(data_find)
data_find=as.data.frame(data_find)
data_find=log(data_find+1)
target=c(rep(0,63),rep(1,77))
data_find=cbind(target,data_find)
data_find$target=as.factor(data_find$target)

data_test=full_join(data23,data25,by="gene")
data_test[is.na(data_test)]<-0
rownames(data_test)=data_test[,1]
data_test=data_test[,-1]
data_test=t(data_test)
data_test=as.data.frame(data_test)
data_test=log(data_test+1)
target=c(rep(0,23),rep(1,25))
data_test=cbind(target,data_test)
data_test$target=as.factor(data_test$target)

data_188[is.na(data_188)]<-0
rownames(data_188)=data_188[,1]
data_188=data_188[,-1]
data_188=t(data_188)
data_188=as.data.frame(data_188)
data_188=log(data_188+1)
target=c(rep(0,23),rep(1,25),rep(0,63),rep(1,77))
data_188=cbind(target,data_188)
data_188[,1]=as.factor(data_188[,1])
#snv位置,只考虑非同义突变##########################################################

for (i in 1:77) {
  find_output_snv_gout[[i]] <- find_output_snv_gout[[i]] %>% filter(mutation_type %in% c("N"))
  find_output_snv_gout[[i]] <- find_output_snv_gout[[i]] %>% select("scaffold","position")
  find_output_snv_gout[[i]]$xmxdhl<-rep(1,nrow(find_output_snv_gout[[i]]))
  find_output_snv_gout[[i]]$tt<- paste(find_output_snv_gout[[i]]$scaffold, find_output_snv_gout[[i]]$position,sep="_")
  find_output_snv_gout[[i]]<-find_output_snv_gout[[i]] %>% select("tt","xmxdhl")
}

for (i in 1:63) {
  find_output_snv_control[[i]] <- find_output_snv_control[[i]] %>% filter(mutation_type %in% c("N"))
  find_output_snv_control[[i]] <- find_output_snv_control[[i]] %>% select("scaffold","position")
  find_output_snv_control[[i]]$xmxdhl<-rep(1,nrow(find_output_snv_control[[i]]))
  find_output_snv_control[[i]]$tt<- paste(find_output_snv_control[[i]]$scaffold,find_output_snv_control[[i]]$position,sep="_")
  find_output_snv_control[[i]]<-find_output_snv_control[[i]] %>% select("tt","xmxdhl")
}


for (i in 1:25) {
  test_output_snv_gout[[i]] <- test_output_snv_gout[[i]] %>% filter(mutation_type %in% c("N"))
  test_output_snv_gout[[i]] <- test_output_snv_gout[[i]] %>% select("scaffold","position")
  test_output_snv_gout[[i]]$xmxdhl<-rep(1,nrow(test_output_snv_gout[[i]]))
  test_output_snv_gout[[i]]$tt<- paste(test_output_snv_gout[[i]]$scaffold,test_output_snv_gout[[i]]$position,sep="_")
  test_output_snv_gout[[i]]<-test_output_snv_gout[[i]] %>% select("tt","xmxdhl")
}

for (i in 1:23) {
  test_output_snv_control[[i]] <-  test_output_snv_control[[i]] %>% filter(mutation_type %in% c("N"))
  test_output_snv_control[[i]] <-  test_output_snv_control[[i]] %>% select("scaffold","position")
  test_output_snv_control[[i]]$xmxdhl<-rep(1,nrow( test_output_snv_control[[i]]))
  test_output_snv_control[[i]]$tt<- paste( test_output_snv_control[[i]]$scaffold, test_output_snv_control[[i]]$position,sep="_")
  test_output_snv_control[[i]]<-test_output_snv_control[[i]] %>% select("tt","xmxdhl")
}
#合并表格##################################################################
for (i in 2:77) {
  find_output_snv_gout[[1]] <-    full_join(find_output_snv_gout[[1]],  find_output_snv_gout[[i]], by = "tt")
}
for (i in 2:63) {
  find_output_snv_control[[1]] <- full_join(find_output_snv_control[[1]], find_output_snv_control[[i]], by = "tt")
}


for (i in 2:25) {
  test_output_snv_gout[[1]] <- full_join(test_output_snv_gout[[1]], test_output_snv_gout[[i]], by = "tt")
}
for (i in 2:23) {
  test_output_snv_control[[1]] <- full_join(test_output_snv_control[[1]], test_output_snv_control[[i]], by = "tt")
}


#去除稀少的snv位置（40%）,减少计算量#############################################################################
snv77= find_output_snv_gout[[1]]
snv63= find_output_snv_control[[1]]
snv25=test_output_snv_gout[[1]]
snv23=test_output_snv_control[[1]]

colnames(snv77)<-c("position",names(find_output_snv_gout))
rownames(snv77)<-snv77[,1]
snv77[is.na(snv77)]<-0

colnames(snv63)<-c("position",names(find_output_snv_control))
rownames(snv63)<-snv63[,1]
snv63[is.na(snv63)]<-0

colnames(snv25)<-c("position",names(test_output_snv_gout))
rownames(snv25)<-snv25[,1]
snv25[is.na(snv25)]<-0

colnames(snv23)<-c("position",names(test_output_snv_control))
rownames(snv23)<-snv23[,1]
snv23[is.na(snv23)]<-0

n1=apply(snv77==1, 1, sum)
i1=which(n1>=30)
snv77i=snv77[i1,]
p1=snv77i$position

n2=apply(snv63==1, 1, sum)
i2=which(n2>=25)
snv63i=snv63[i2,]
p2=snv63i$position

n3=apply(snv25==1, 1, sum)
i3=which(n3>=10)
snv25i=snv25[i3,]
p3=snv25i$position

n4=apply(snv23==1, 1, sum)
i4=which(n4>=9)
snv23i=snv23[i4,]
p4=snv23i$position

pp=union(p1,p2)
pp2=union(p3,p4)
pp3=union(pp,pp2)
#整合表格##################################################

snv140=full_join(snv63,snv77,by="position")
snv48=full_join(snv23,snv25,by="position")
snv188=full_join(snv140,snv48,by="position")

snv188=subset(snv188,position%in%pp3)
rownames(snv188)=snv188[,1]
snv188=snv188[,-1]
snv188[is.na(snv188)]<-0
snv188t=t(snv188)
snv188t=as.data.frame(snv188t)
target_snv=c(rep(0,63),rep(1,77),rep(0,23),rep(1,25))
snv188t=cbind(target_snv,snv188t)
snv188t[,1:ncol(snv188t)]=lapply(snv188t[,1:ncol(snv188t)], factor)

snv140=subset(snv140,position%in%pp3)
rownames(snv140)=snv140[,1]
snv140=snv140[,-1]
snv140[is.na(snv140)]<-0
snv140t=t(snv140)
snv140t=as.data.frame(snv140t)
target_snv=c(rep(0,63),rep(1,77))
snv140t=cbind(target_snv,snv140t)
snv140t[,1:ncol(snv140t)]=lapply(snv140t[,1:ncol(snv140t)], factor)


snv48=subset(snv48,position%in%pp3)
rownames(snv48)=snv48[,1]
snv48=snv48[,-1]
snv48[is.na(snv48)]<-0
snv48t=t(snv48)
snv48t=as.data.frame(snv48t)
target_snv=c(rep(0,23),rep(1,25))
snv48t=cbind(target_snv,snv48t)
snv48t[,1:ncol(snv48t)]=lapply(snv48t[,1:ncol(snv48t)], factor)























