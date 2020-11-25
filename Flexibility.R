#Lecane inermis & Daphnia magna evaluation of flexibility of microbiome

library("vegan")
library("RVAideMemoire")
Data <- read.csv("otu_table.csv") 
row.names(Data)<-Data[,1]
Data<-Data[,-1]
OTUList<-as.data.frame(Data)
rm(Data)
tOTUList<-as.data.frame(t(OTUList))
taxaall<-read.csv("taxonomy.csv")
taxaNC<-subset(taxaall, taxaall$Class!="c:Chloroplast")
taxaNM<-subset(taxaall, taxaall$Family!="f:Mitochondria")
taxa<-subset(taxaNC, taxaNC$Phylum!="")
ttaxa<-as.data.frame(t(taxa))
rownames(taxa)<-taxa[,1]
OTUnC<-subset(OTUList, taxaall$Class!="c:Chloroplast")
OTU<-subset(OTUnC, taxaNC$Phylum!="")
tOTU<-t(OTU)
reads<-as.data.frame(rowSums(tOTU))
tOTU<-as.data.frame(t(OTU))
tOTUList<-subset(tOTU, reads[,1]>2000)
# variables<-subset(variables, variables$X.OTU.ID%in%row.names(tOTUList))
OTUList<-as.data.frame(t(tOTUList))
variExp<-read.csv("vari_exp.csv")
tOTUList<-as.data.frame(t(OTUList))
Exp<-subset(tOTUList, row.names(tOTUList)%in%variExp$Name)
variExp<-subset(variExp,variExp$Name%in%row.names(tOTUList))


tExp<-as.data.frame(t(Exp))
reads<-rowSums(Exp)


library("GUniFrac")
rareOTU<-Rarefy(Exp, depth = min(reads))  #Rarefaction
min(reads)
rffOTU<-as.data.frame(rareOTU$otu.tab.rff)  #5435
totalSeq<-rowSums(rffOTU)
trffOTU<-t(rffOTU)
raOTU <- rffOTU[ ,colSums(rffOTU)!=0]
traOTU<-t(raOTU)
paOTU<-raOTU
paOTU[paOTU>0]=1
row.names(raOTU)<-variExp$Short

beta<-vegdist(raOTU, method="bray")
cave <- hclust(beta, method="average")
plot(cave, main="bray curtis", hang=-1, cex=0.7)
pairwise.perm.manova(beta, variExp$aniwat)

rot<-subset(raOTU, variExp$GENOT=="rot")
variRot<-subset(variExp, variExp$GENOT=="rot")
rotR<-subset(rot, variRot$aniwat=="rotifer")
Vrot<-subset(variRot,variRot$aniwat=="rotifer")

beta<-vegdist(rotR, method="bray")
cave <- hclust(beta, method="complete") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)

adonis(beta~Vrot$species*Vrot$Food)

 
LHT1<-subset(rotR, Vrot$Exp_temp=="Y")
LHT<-LHT1[ ,colSums(LHT1)>0]
VLHT<-subset(Vrot, Vrot$Exp_temp=="Y")
betaT<-vegdist(LHT, method="bray")
cave <- hclust(betaT, method="average") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)
adonis(betaT~VLHT$Temp+VLHT$Food+VLHT$Original_Food)

LHTi1<-subset(rotR, Vrot$Exp_time=="Y")
LHTi<-LHTi1[ ,colSums(LHTi1)>0]
VLHTi<-subset(Vrot, Vrot$Exp_time=="Y")
betaT<-vegdist(LHTi, method="bray")
cave <- hclust(betaT, method="complete") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)
adonis(betaT~VLHTi$Time+VLHTi$Food)


daph<-subset(raOTU, variExp$species=="Daphnia magna")
varidaph<-subset(variExp, variExp$species=="Daphnia magna")
beta<-vegdist(daph, method="bray")
cave <- hclust(beta, method="complete") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)

adonis(beta~varidaph$place)


daphT<-subset(daph, varidaph$Exp_time=="Y")
varidaphT<-subset(varidaph, varidaph$Exp_time=="Y")

beta<-vegdist(daphT, method="bray")
cave <- hclust(beta, method="complete") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)

adonis(beta~varidaphT$GENOT+varidaphT$Time)

daphTemp<-subset(daph, varidaph$Exp_temp=="Y")
varidaphTemp<-subset(varidaph, varidaph$Exp_temp=="Y")

beta<-vegdist(daphTemp, method="bray")
cave <- hclust(beta, method="complete")
plot(cave, main="bray curtis", hang=-1, cex=0.7)

adonis(beta~varidaphTemp$Temp)

library("ggplot2")
library("reshape2")
#plot(beta~varidaph)
betam<-as.data.frame(as.matrix(beta))
#betamm<-melt(betam)

groups <- varidaph$place
library("betapart")

par(mfrow=c(1,2))

dist<-bray.part(daph)
bd<-betadisper(dist[[3]], groups)
df <- data.frame(Distance_to_centroid=bd$distances,Group=bd$group)
plot(Distance_to_centroid~Group, data=df)

groups <- variRot$species
dist<-bray.part(rot)
bd<-betadisper(dist[[3]], groups)
df <- data.frame(Distance_to_centroid=bd$distances,Group=bd$group)
plot(Distance_to_centroid~Group, data=df)





ggplot(data=df,aes(x=Group,y=Distance_to_centroid,colour=groups)) +
  geom_boxplot(alpha=0.5) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))


####the same for water



rotW<-subset(rot, variRot$aniwat=="water")
VWat<-subset(variRot,variRot$aniwat=="water")

beta<-vegdist(rotW, method="bray")
cave <- hclust(beta, method="complete") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)

adonis(beta~VWat$species*VWat$Food)


LHT1<-subset(rotW, VWat$Exp_temp=="Y")
LHT<-LHT1[ ,colSums(LHT1)>0]
VLHT<-subset(VWat, VWat$Exp_temp=="Y")
betaT<-vegdist(LHT, method="bray")
cave <- hclust(betaT, method="complete") 
plot(cave, main="bray curtis", hang=-1, cex=0.7)
adonis(betaT~VLHT$Temp+VLHT$Food+VLHT$Original_Food)

LHTi1<-subset(rotW, VWat$Exp_time=="Y")
LHTi<-LHTi1[ ,colSums(LHTi1)>0]
VLHTi<-subset(VWat, VWat$Exp_time=="Y")
betaT<-vegdist(LHTi, method="bray")
cave <- hclust(betaT, method="complete")
plot(cave, main="bray curtis", hang=-1, cex=0.7)
adonis(betaT~VLHTi$Time+VLHTi$Food)





