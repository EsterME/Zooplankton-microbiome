#Lecane inermis & Daphnia magna evaluation of flexibility of microbiome
setwd("C:/Users/Admin/Desktop/CNR/Rotifers/Results_Sept2019/")
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

tOTU1<-as.data.frame(t(OTU))
vari<-read.csv("variables_new_nmds.csv")
vari2<-subset(vari, vari$X%in%row.names(tOTU1))
tOTU<-subset(tOTU1, row.names(tOTU1)%in%vari2$X)
tOTUListA<-subset(tOTU, vari2$substrate=="animal")
variA<-subset(vari2, vari2$substrate=="animal")
reads<-as.data.frame(rowSums(tOTUListA))
tOTUhigh<-subset(tOTUListA, reads[,1]>2000)
tOTUlow<-subset(tOTUListA, reads[,1]<2000)
varih<-subset(variA, variA$X%in%row.names(tOTUhigh))
varil<-subset(variA, variA$X%in%row.names(tOTUlow))


OTUhigh<-as.data.frame(t(tOTUhigh))
readsH<-rowSums(tOTUhigh)

library("GUniFrac")
rareOTU<-Rarefy(tOTUListA, depth = min(readsH))  #Rarefaction
min(readsH)
rffOTU<-as.data.frame(rareOTU$otu.tab.rff)  #2207
totalSeq<-rowSums(rffOTU)
trffOTU<-t(rffOTU)
raOTU<-rbind(rffOTU,tOTUlow)
raOTU <- raOTU[ ,colSums(rffOTU)!=0]
variables<-rbind(varih, varil)
row.names(raOTU)<-variables$uni
traOTU<-t(raOTU)
paOTU<-raOTU
paOTU[paOTU>0]=1


beta<-vegdist(raOTU, method="bray")
cave <- hclust(beta, method="complete")
plot(cave, main="bray curtis", hang=-1, cex=0.7)
adonis(beta~variables$species)

library("RColorBrewer")
NMDS<-metaMDS(raOTU, distance="bray", k=3)
plot(NMDS)
col8=brewer.pal(12, "Paired")
colNew<-c('#bcf60c',	'#3cb44b',	'#808000',		"thistle3",	'#f58231',	'#9a6324',	"tan1",	"slategray4",	'grey78',	'#000000',	"snow4",	"slategray2",	'#fffac8',	'#f032e6',	"mediumorchid",	'#800000','#e6194b',	'#ffe119',	"pink3",	'#fabebe',	'#ffd8b1',"peachpuff3",	'#aaffc3',	'steelblue1',	'#4363d8',	'#46f0f0',	'#e6beff',	'#911eb4')
shapeNew<-c("22","7","24","23","21","10")
data.scores <- as.data.frame(scores(NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.factor(variables$spec_place)  #  add the grp variable created earlier
data.scores$taxa <- as.factor(variables$shape)
#head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(NMDS, "sites"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=grp,fill=grp, shape=taxa),size=5) + # add the point markers
  scale_fill_manual(values=colNew) +
  scale_color_manual(values=colNew)+
  scale_shape_manual(values=c(23,21,24,22))+
  coord_equal() +
  theme_bw()
#



