library("RAM")
setwd("C:/Users/Admin/Desktop/CNR/Rotifers/Results_Sept2019/")
Data <- read.OTU("otu_table_RAM.csv", sep=",") 
otu <- read.csv("otu_table.csv") 
row.names(otu)<-otu[,1]
otu<-otu[,-1]
variables<-read.meta("variables_new2.csv", sep=",")
tData<-t(Data)
tOTUList<-subset(tData, row.names(tData)%in%row.names(variables))
variAB<-subset(variables, row.names(variables)%in%row.names(tOTUList))
tOTUind<-as.data.frame(subset(tOTUList, variAB$indval!="NA"))
variInd<-as.data.frame(subset(variAB, variAB$indval!="NA"))
OTUind<-as.data.frame(t(tOTUind))
taxonomy<-Data$taxonomy
OTUind2<-cbind(OTUind,taxonomy)
taxa<-read.csv("taxonomyRAM.csv")
rownames(taxa)<-taxa[,1]
taxa<-taxa[,-1]
variables=variInd


####DATA EXPERIMENT######

library("GUniFrac")
#


Species<-c("Adineta vaga", "Lecane inermis", "Rotaria macrura", "Daphnia magna")
vabSp<-variables[variables$species%in%Species,]
dataSp<-Data[,colnames(Data)%in%rownames(vabSp)]
dataSpt<-cbind(dataSp, Data$taxonomy)
colnames(dataSpt)[62]<-"taxonomy"

data=list(data=dataSpt)
input.ranks <- c("kingdom", "phylum", "class", "order",
                 "family", "genus")
ndata <- reformat.taxonomy(data,input.ranks=input.ranks,
                           sep=",")

corRot<-core.OTU(ndata, vabSp, meta.factor="species", percent=0.666665)




############
common<-Reduce(intersect, list(corRot$data$`Adineta vaga`$otuID, corRot$data$water$otuID))
nt<-Data[,-88]

tnt<-t(nt)



##############
#Unifrac dist
#############


# #Reduce replicates
# # First find the maxvalues
Data <- read.csv("otu_table2.csv") 
row.names(Data)<-Data[,1]
Data<-Data[,-1]
OTUList<-as.data.frame(Data)
#rm(Data)
tOTUList<-as.data.frame(t(OTUList))
variables<-read.csv("variables_new2.csv")
tOTUList<-subset(tOTUList, row.names(tOTUList)%in%variables$X)
OTUList<-t(tOTUList)
variables<-subset(variables, variables$X%in%row.names(tOTUList))
taxaall<-read.csv("taxonomy2.csv")
taxaNC<-subset(taxaall, taxaall$Class!="c:Chloroplast")
taxaNM<-subset(taxaall, taxaall$Family!="Mitochondria")
taxa<-subset(taxaNC, taxaNC$Phylum!="")
ttaxa<-as.data.frame(t(taxa))
rownames(taxa)<-taxa[,1]
taxa<-taxa[,-1]
OTUnC<-subset(OTUList, taxaall$Class!="c:Chloroplast")
OTU<-as.data.frame(subset(OTUnC, taxaNC$Phylum!=""))
OTU3<-dataSp
tOTU3<-t(OTU3)


reads<-colSums(OTU3)
vab4<-vabSp
vab4$reads<-reads
List_max<-as.data.frame(tapply(vab4$reads,vab4$species,  max))
variables_red<-vab4$reads%in%(tapply(vab4$reads,vab4$species,  max))
variNR<-subset(vab4, variables_red==TRUE)
# 
OTUNR<-as.data.frame(tOTU3[row.names(tOTU3)%in%row.names(variNR),] )

rownames(OTUNR)<-variNR$species
tOTUNR<-t(OTUNR)
specRep<-c("Adineta vaga", "Daphnia magna",   "Lecane inermis" , "Rotaria macrura" )
tOTUo <- tOTUNR[,specRep]
paOTU<-tOTUo
paOTU[paOTU>0]=0
OTUc<-data.frame()
library("dplyr")
for (i in 1:4) {
  
  coreOTU<-as.data.frame(corRot$data[[i]]$otuID)
  paOTU[,i][rownames(paOTU)%in%coreOTU$`corRot$data[[i]]$otuID`]=1
  OTUc<-rbind(OTUc, coreOTU)
}
list_otu<-unique(OTUc$`corRot$data[[i]]$otuID`)

library(ape)
treeOTU<-read.nexus("raxml_rootedotu.nexus")


OTU2<-paOTU[rownames(paOTU)%in%list_otu,]
tOTU2<-t(OTU2)
tips<-colnames(tOTU2)
tree2<-keep.tip(treeOTU,tips)
Unidist<-GUniFrac(tOTU2, tree2)$unifracs
dw <- Unidist[, , "d_UW"] # Weighted UniFra
cldw<-hclust(as.dist(dw), method="average")
plot(cldw, hang=-1)

OTU2<-OTU2[,cldw$order]

Species<-c("Adineta vaga" ,  "Daphnia magna", "Lecane inermis" , "Rotaria macrura"   )
fam<-as.data.frame(sort(unique(taxa$Family)))
colnames(fam)<-"X1"
for (i in 1:4) { 
  spec<-Species[[i]]
  ntLE<-as.data.frame(subset(tOTU3,vab4$species==spec))
  PAntLE<-ntLE
  PAntLE[PAntLE>0]=1
  richLE<-as.data.frame(rowSums(PAntLE))
  corLE<-length(corRot$data[[i]]$otuID)
  
  #############################Loop for all species
  tntLE<-as.data.frame(t(ntLE))
  perLE<-NULL
  #######################################
  for (j in colnames(tntLE) ){
    OTU<-tntLE[[j]]!=0
    OTUL<-as.list(subset(rownames(tntLE), OTU==TRUE))
    common<-Reduce(intersect, list(corRot$data[[i]]$otuID, OTUL))
    perc<-length(common)/length(OTUL)
    perLE <- as.data.frame(rbind(perLE, perc))
  }
  perLE[,2]<-spec
  print(perc)
  
  write.table(perLE, "percentageSpec.txt", append = T)
  ########################################
  spec<-Species[[i]]
  otus<-corRot$data[[i]]$otuID
  taxCore<-subset(taxa, rownames(taxa)%in%otus)
  taxCore[,7]<-rep(1, length(taxCore[,1]))
  sumTaxa1<-as.list(by(taxCore$V7,taxCore$Family, sum))
  Tax1<-as.data.frame(unlist(sumTaxa1))
  Tax1$X1<-row.names(Tax1)
  fam<-left_join(fam,Tax1,by="X1")
  colnames(fam)[i+1]<-spec
}



library("reshape2")
library("ggplot2")
rownames(fam)<-fam[,1]
fam<-fam[,-1]
fam[is.na(fam)] <- 0
fam2<-subset(fam, rowSums(fam)!=0)

#library("textshape")
colS<-colSums(fam2)
famR<-t(fam2)/colSums(fam2)
tfamR<-as.data.frame(t(famR*100))
famM<-as.matrix(famR)
tfamRo <- tfamR[,cldw$order]
mfam<-melt(cbind(tfamRo, ind=rownames(tfamRo), id.vars=c('ind')) )

ggplot(mfam,aes(x = variable, y = ind ,size = value)) + geom_point(alpha=0.4) +scale_size(limits = c(0.1,100)) +theme_light()

#######################################################################
###################Percentage core microbiome##########################
#######################################################################
#######

spec<-"Rotaria macrura"
ntRM<-subset(tnt,variables$species==spec)

PAntRM<-ntRM
PAntRM[PAntRM>0]=1
richRM<-as.data.frame(rowSums(PAntRM))
corRM<-length(corRot$data$`Rotaria macrura`$otuID)

#############################Loop for all species
tntRM<-as.data.frame(t(ntRM))
perRM<-NULL
for (i in colnames(tntRM) ){
  OTU<-tntRM[[i]]!=0
  OTUL<-as.list(subset(rownames(tntRM), OTU==TRUE))
  common<-Reduce(intersect, list(corRot$data$`Rotaria macrura`$otuID, OTUL))
  perc<-length(common)/length(OTUL)
  perRM <- as.data.frame(rbind(perRM, perc))
  print(perc)
  
}

perRM[,2]<-spec
colnames(perRM)<-c("percentage", "species")
###########################

spec<-"Lecane inermis"
ntLE<-subset(tnt,variables$species==spec)
PAntLE<-ntLE
PAntLE[PAntLE>0]=1
richLE<-as.data.frame(rowSums(PAntLE))
corLE<-length(corRot$data$`Lecane inermis`$otuID)

#############################Loop for all species
tntLE<-as.data.frame(t(ntLE))
perLE<-NULL
for (i in colnames(tntLE) ){
  OTU<-tntLE[[i]]!=0
  OTUL<-as.list(subset(rownames(tntLE), OTU==TRUE))
  common<-Reduce(intersect, list(corRot$data$`Lecane inermis`$otuID, OTUL))
  perc<-length(common)/length(OTUL)
  perLE <- as.data.frame(rbind(perLE, perc))
  print(perc)
  
}

perLE[,2]<-spec
colnames(perLE)<-c("percentage", "species")
###########################


spec<-"Daphnia magna"
ntDa<-subset(tnt,variables$species==spec)
PAntDa<-ntDa
PAntDa[PAntDa>0]=1
richDa<-as.data.frame(rowSums(PAntDa))
corDa<-length(corRot$data$`Daphnia magna`$otuID)

#############################Loop for all species
tntDa<-as.data.frame(t(ntDa))
perDa<-NULL
for (i in colnames(tntDa) ){
  OTU<-tntDa[[i]]!=0
  OTUL<-as.list(subset(rownames(tntDa), OTU==TRUE))
  common<-Reduce(intersect, list(corRot$data$`Daphnia magna`$otuID, OTUL))
  perc<-length(common)/length(OTUL)
  perDa <- as.data.frame(rbind(perDa, perc))
  print(perc)
  
}

perDa[,2]<-spec
colnames(perDa)<-c("percentage", "species")






###########################
hammer<-read.csv("hammer_data.csv")
perLE<-read.csv("percentageSpec.csv")
perc<-rbind(perLE, hammer)
x<-perc$species
print(levels(x))  ## This will show the levels of x are "Levels: a b c d e"
spec<-c("Daphnia magna", "Adineta vaga", "Lecane inermis" , "Rotaria macrura",  "caterpillars", "other animals"  )

perc$species <- factor(perc$species , levels=c("Daphnia magna", "Adineta vaga", "Lecane inermis" , "Rotaria macrura",  "caterpillars", "other animals"  ))
boxplot(perc$percentage~perc$species, main="species level", ylab="% core of reads", xlab="species")



###################################################
#Genus core:
Data2 <- read.OTU("otu_table_RAM.csv", sep=",") 
Data2<-subset(Data2, rownames(Data2)%in%rownames(taxa))
Data<-read.csv("otu_table2.csv")
rownames(Data)<-Data[,1]
Data<-Data[,-1]
Data<-subset(Data, rownames(Data)%in%rownames(taxa))
tData<-as.data.frame(t(Data))
tData<-as.data.frame(subset(tData, rownames(tData)%in%variables$X))
tDataN<-as.data.frame(subset(tData, variables$substrate=="animal"))
vari2<-subset(variables, variables$X%in%rownames(tDataN) )
DataN<-as.data.frame(t(tDataN))
###Loop for replicates
md4<-as.data.frame((summary(as.factor(vari2$spec_place)))>=3) #find species with more than 3 replicates
colnames(md4)<-"count"
spec4<-subset(md4, md4$count==TRUE)
vab4<-subset(vari2, vari2$spec_place%in%rownames(spec4))
data4r<-DataN[,colnames(DataN)%in%vab4$X]
data4<-cbind(data4r, Data2$taxonomy)
rownames(vab4)<-vab4$X
colnames(data4)[71]<-"taxonomy"
data=list(data=data4)

input.ranks <- c("kingdom", "phylum", "class", "order", 
                 "family", "genus")
ndata <- reformat.taxonomy(data,
                           input.ranks=input.ranks,
                           sep=",")
valid.taxonomy(ndata)
coreGen<-core.OTU(ndata, vab4, meta.factor="spec_place",percent = 0.6667)

tdat4<-as.data.frame(t(data4r))
specRep<-c("CopepodLM","Daphnia galeataLM", "Daphnia magnadecaeste", "Daphnia magnaebert", "Daphnia magnaFinnland", "Daphnia obtusIRSA", "large calanoid copepodA993", "Lecane inermisLIH", "Lecane inermisLIK", "Rotaria macruraSt14")
fam<-as.data.frame(sort(unique(taxa$Family)))

OTU3<-data4r
tOTU3<-t(OTU3)


reads<-colSums(OTU3)

vab4$reads<-reads
List_max<-as.data.frame(tapply(vab4$reads,vab4$spec_place,  max))
variables_red<-vab4$reads%in%(tapply(vab4$reads,vab4$spec_place,  max))
variNR<-subset(vab4, variables_red==TRUE)
# 
OTUNR<-as.data.frame(tOTU3[row.names(tOTU3)%in%rownames(variNR),] )
rownames(OTUNR)<-variNR$spec_place
OTU<-OTUNR[-3,]
tOTU<-as.data.frame(t(OTU))
specRep<-c("CopepodLM","Daphnia galeataLM", "Daphnia magnadecaeste", "Daphnia magnaebert", "Daphnia magnaFinnland", "Daphnia obtusIRSA", "large calanoid copepodA993", "Lecane inermisLIH", "Lecane inermisLIK", "Rotaria macruraSt14")

tOTUo <- tOTU[,specRep]
paOTU<-tOTUo
paOTU[paOTU>0]=0
OTUc<-data.frame()
for (i in 1:10) {
  coreOTU<-as.data.frame(coreGen$data[[i]]$otuID)
  paOTU[,i][rownames(paOTU)%in%coreOTU$`coreGen$data[[i]]$otuID`]=1
  OTUc<-rbind(OTUc, coreOTU)
}
list_otu<-unique(OTUc$`coreGen$data[[i]]$otuID`)

library(ape)
treeOTU<-read.nexus("raxml_rootedotu.nexus")


OTU2<-paOTU[rownames(paOTU)%in%list_otu,]
OTU2<-paOTU[rowSums(paOTU)!=0,]
tOTU2<-t(OTU2)
tips<-rownames(OTU2)
tree2<-keep.tip(treeOTU,tips)
Unidist<-GUniFrac(tOTU2, tree2)$unifracs
dw <- Unidist[, , "d_UW"] # Unweighted UniFra
cldw<-hclust(as.dist(dw), method="average")
plot(cldw, hang=-1)

OTU3<-OTU2[,cldw$order]
fam<-as.data.frame(sort(unique(taxa$Family)))
colnames(fam)<-"X1"

for (i in 1:10) { 
  spec<-specRep[[i]]
  ntLE<-subset(tdat4,vab4$spec_place==spec)
  PAntLE<-ntLE
  PAntLE[PAntLE>0]=1
  richLE<-as.data.frame(rowSums(PAntLE))
  corLE<-length(coreGen$data[[i]]$otuID)
  
  #############################Loop for all species
  tntLE<-as.data.frame(t(ntLE))
  perLE<-NULL
  #######################################
  for (j in colnames(tntLE) ){
    OTU<-tntLE[[j]]!=0
    OTUL<-as.list(subset(rownames(tntLE), OTU==TRUE))
    common<-Reduce(intersect, list(coreGen$data[[i]]$otuID, OTUL))
    perc<-length(common)/length(OTUL)
    perLE <- as.data.frame(rbind(perLE, perc))
  }
  perLE[,2]<-spec
  print(perc)
  
  #write.table(perLE, "percentageRep.txt", append = T)
  
  library("reshape2")
  library("ggplot2")
  ########################################
  
  spec<-specRep[[i]]
  otus<-coreGen$data[[i]]$otuID
  taxCore<-subset(taxa, rownames(taxa)%in%otus)
  taxCore[,7]<-rep(1, length(taxCore[,1]))
  sumTaxa1<-as.list(by(taxCore$V7,taxCore$Family, sum))
  Tax1<-as.data.frame(unlist(sumTaxa1))
  Tax1$X1<-row.names(Tax1)
  fam<-left_join(fam,Tax1,by="X1")
  colnames(fam)[i+1]<-spec
  
}

rownames(fam)<-fam[,1]
fam<-fam[,-1]
fam[is.na(fam)] <- 0
fam2<-subset(fam, rowSums(fam)!=0)
colS<-colSums(fam2)
famR<-t(fam2)/colSums(fam2)
tfamR<-as.data.frame(t(famR*100))
famM<-as.matrix(famR)
specRepS<-c("Lecane inermisLIH", "Lecane inermisLIK", "Rotaria macruraSt14","large calanoid copepodA993", "CopepodLM","Daphnia galeataLM", "Daphnia magnadecaeste", "Daphnia magnaebert", "Daphnia magnaFinnland", "Daphnia obtusIRSA")
tfamRo <- tfamR[,cldw$order]
mfam<-melt(cbind(tfamRo, ind=rownames(tfamRo), id.vars=c('ind')) )

ggplot(mfam,aes(x = variable, y = ind ,size = value)) + geom_point(alpha=0.4) +scale_size(limits = c(0.1,100)) +theme_light()

rep<-read.csv("percentageRepread.csv")
spec<-c("Adineta vaga", "Lecane inermis" , "Rotaria macrura", "Daphnia magna", "caterpillars", "other animals"  )
spec<-unique(as.character(rep$species))
rep$species <- factor(rep$species , levels=c("R. mac","E. pada",  "D. gal",  "D. mag 1", "D. mag 2", "D. mag 3", "D. obt ", "cal cop" ,"L. in H" , "L. in K"))

plot(rep$percentage~rep$species, main="replicates")




