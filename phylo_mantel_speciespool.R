setwd("C:/Users/Admin/Desktop/CNR/Rotifers/Results_Sept2019/")
library("vegan")
library("GUniFrac")
# #Reduce replicates
# # First find the maxvalues
Data <- read.csv("otu_table2.csv") 
row.names(Data)<-Data[,1]
Data<-Data[,-1]
OTUList<-as.data.frame(Data)
rm(Data)
tOTUList<-as.data.frame(t(OTUList))
variables<-read.csv("variables_new2.csv")
tOTUList<-subset(tOTUList, row.names(tOTUList)%in%variables$X)
OTUList<-t(tOTUList)
variables<-subset(variables, variables$X%in%row.names(tOTUList))
taxaall<-read.csv("taxonomy.csv")
taxaNC<-subset(taxaall, taxaall$Class!="c:Chloroplast")
taxaNM<-subset(taxaall, taxaall$Family!="f:Mitochondria")
taxa<-subset(taxaNC, taxaNC$Phylum!="")
ttaxa<-as.data.frame(t(taxa))
rownames(taxa)<-taxa[,1]
OTUnC<-subset(OTUList, taxaall$Class!="c:Chloroplast")
OTU<-subset(OTUnC, taxaNC$Phylum!="")
tOTU<-t(OTU)
tOTU <- as.data.frame(tOTU[ ,colSums(tOTU)!=0] )
OTU<-t(tOTU)
species<-list(variables$species)
aggOTU1<-with(tOTU, aggregate(tOTU, by=species, FUN=sum))
EpiLecE<-aggOTU1[c(9,15),-1]
rownames(EpiLecE)<-c("Epiphanes_senta_A962ES"  ,  "Lecane_elsa_A984LE" )
library("GUniFrac")
rareOTU<-Rarefy(tOTU, depth = 1000)  
raOTU<-as.data.frame(rareOTU$otu.tab.rff)  #5435

traOTU<-as.data.frame(t(raOTU))
variables2<-subset(variables, variables$X%in%row.names(raOTU))
phyOTU<-as.data.frame(traOTU)
tphyOTU<-as.data.frame(t(phyOTU))
species<-list(variables2$species)
aggOTU<-with(tphyOTU, aggregate(tphyOTU, by=species, FUN=sum))
phyloname<-c("Adineta_vaga_culture" , "Brachionus_quadridentatus_A988BQ", "cyclopid_A984C","Daphnia_galeata_DLM","Daphnia_magna","Daphnia_obtusa","Diaphanosoma_brachyurum_A995C", "Euchlanis_dilatata_A991E","Eudiaptomus_padanus_CLM","Keratella_serrulata_A990K","Keratella_quadrata_A985K", "calanoid_large_A993C", "Lecane_inermis_culture","Mesocyclops_leuckarti_A995CO", "calanoid_nauplii_A986C","Polyarthra_dolichoptera_A989P","Rotaria_macrura_A970RM","Rotaria_rotatoria_A994RR","Simocephalus_sp_A994S","Simocephalus_vetulus_A992SP","water" )
rownames(aggOTU)<-phyloname
aggOTU<-aggOTU[,-1]
rareOTU<-Rarefy(aggOTU, depth = min(rowSums(aggOTU)))  
min(rowSums(aggOTU))#Rarefaction
rffOTU<-as.data.frame(rareOTU$otu.tab.rff)  #5435
rffOTU<-rbind(rffOTU, EpiLecE)
raOTU <- rffOTU[ ,colSums(rffOTU)!=0]
traOTU<-as.data.frame(t(raOTU))
treeZoo<-read.nexus("BEAST_phylo_new.nex")
vari_phyl<-read.csv("vaPhy.csv")
tipsL<-treeZoo$tip.label
phyOTU2<-raOTU[rownames(raOTU)%in%tipsL,]
phyOTU2<-phyOTU2[order(row.names(phyOTU2)),]
variphy<-read.csv("variPHY.csv")
#treeZoo2<-drop.tip(treeZoo, tips) #, trim.internal = TRUE, subtree = FALSE,
#root.edge = 0, rooted = is.rooted(treeZoo), collapse.singles = TRUE,
#interactive = FALSE)
pruned.tree<-drop.tip(treeZoo,treeZoo$tip.label[-match(tipsL, treeZoo$tip.label)])
dist_treeDi<-as.dist(cophenetic.phylo(pruned.tree))
library("vegan")
dist_beta<-vegdist(phyOTU2, method="bray", diag=TRUE, binary = F)
cave <- hclust(dist_beta, method="average")
plot(cave, main="bray curtis", hang=-1, cex=0.7)

mantel(dist_treeDi, dist_beta)

dist_beta<-vegdist(phyOTU2, method="jaccard", diag=TRUE, binary=TRUE)
cave <- hclust(dist_beta, method="average")
plot(cave, main="bray curtis", hang=-1, cex=0.7)

mantel(dist_treeDi, dist_beta)

library("betapart")
paphy<-phyOTU2
paphy[paphy>0]=1
dist_sor<-beta.pair(paphy)

mantel(dist_treeDi, dist_sor$beta.sor)

vari_phyl<-read.csv("vaPhy.csv")
vari_phyl<-vari_phyl[order(vari_phyl$Phylo_name),]
vari_phyl_Cr<-subset(vari_phyl, vari_phyl$rot_crus=="crustacean")
tipsLC<-as.vector(vari_phyl_Cr$Phylo_name)
pruned.treeC<-drop.tip(treeZoo,treeZoo$tip.label[-match(tipsLC, treeZoo$tip.label)])
dist_treeDiC<-as.dist(cophenetic.phylo(pruned.treeC))

phyC<-as.data.frame(subset(phyOTU2, vari_phyl$rot_crus=="crustacean"))
distC<-vegdist(phyC, method="bray", diag=TRUE)
cave <- hclust(distC, method="average") #average linkage !Seems the best one for our data!
plot(cave, main="bray curtis", hang=-1, cex=0.7)
mantel(dist_treeDiC, distC)
distC<-vegdist(phyC, method="bray", diag=TRUE, binary=T)
mantel(dist_treeDiC, distC)

vari_phyl_Cr<-subset(vari_phyl, vari_phyl$rot_crus=="rotifer")
tipsLC<-as.vector(vari_phyl_Cr$Phylo_name)
pruned.treeC<-drop.tip(treeZoo,treeZoo$tip.label[-match(tipsLC, treeZoo$tip.label)])
dist_treeDiC<-as.dist(cophenetic.phylo(pruned.treeC))
phyC<-as.data.frame(subset(phyOTU2, vari_phyl$rot_crus=="rotifer"))
distC<-vegdist(phyC, method="bray", diag=TRUE)
cave <- hclust(distC, method="average") #average linkage !Seems the best one for our data!
plot(cave, main="bray curtis", hang=-1, cex=0.7)
mantel(distC, dist_treeDiC)
distC<-vegdist(phyC, method="bray", diag=TRUE, binary=T)
mantel(distC, dist_treeDiC)

vari_phyl_Cr<-subset(vari_phyl2, vari_phyl2$Taxa=="cop")
tipsLC<-as.vector(vari_phyl_Cr$Phylo_name)
pruned.treeC<-drop.tip(treeZoo,treeZoo$tip.label[-match(tipsLC, treeZoo$tip.label)])
dist_treeDiC<-as.dist(cophenetic.phylo(pruned.treeC))

phyC<-as.data.frame(subset(phyOTU2, vari_phyl2$Taxa=="cop"))
distC<-vegdist(phyC, method="bray", diag=TRUE)
cave <- hclust(distC, method="average") #average linkage !Seems the best one for our data!
plot(cave, main="bray curtis", hang=-1, cex=0.7)
mantel(distC, dist_treeDiC)


betaO<-vegdist(overRep, method="jaccard")


mantel(dist_beta,dist_treeDi)
mantel(betaO,dist_treeDi)

library("GUniFrac")
library("ape")
ptreeZotu<-read.nexus("tree_new.nexus")
phyUni<-read.csv("phyOTU2.csv")
row.names(phyUni)<-phyUni[,1]
phyUni<-phyUni[,-1]
tphyUni<-as.data.frame(t(phyUni))

lOTU<-as.list(rownames(tphyUni))
pruned.tree<-drop.tip(ptreeZotu,ptreeZotu$tip.label[-match(lOTU, ptreeZotu$tip.label)])
Uni<-GUniFrac(phyUni, pruned.tree)
unif<-Uni$unifracs
dw <- as.dist(unif[, , "d_VAW"])
cave <- hclust(dist_beta, method="average") #average linkage !Seems the best one for our data!
plot(cave, main="abund weighted Unifrac", hang=-1, cex=0.7)

mantel(dw,dist_treeDi)

dw <- as.dist(unif[, , "d_UW"])
cave <- hclust(dist_beta, method="average") #average linkage !Seems the best one for our data!
plot(cave, main="abund weighted Unifrac", hang=-1, cex=0.7)

mantel(dw,dist_treeDi)
