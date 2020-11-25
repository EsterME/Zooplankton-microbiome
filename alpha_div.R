setwd("C:/Users/Admin/Desktop/CNR/Rotifers/Results_Sept2019/")


Data <- read.csv("otu_table.csv") 
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
reads<-as.data.frame(rowSums(tOTU))
tOTU<-as.data.frame(t(OTU))
tOTUList<-subset(tOTU, reads[,1]>2000)

library("ape")
list_otu<-row.names(OTU)
fasta_zotu<-read.FASTA("zotus.fa", type="DNA")
fasta_taxa<-subset(fasta_zotu, fasta_zotu%in%list_otu)
fasta_taxa<-fasta_zotu[names(fasta_zotu) %in% list_otu]
write.FASTA(fasta_taxa, "fasta_indentified.fasta", header = NULL, append = FALSE)

variables<-subset(variables, variables$X%in%row.names(tOTUList))
OTUList<-as.data.frame(t(tOTUList))
#POSTHOC FOR ADONIS:
#Bray-Curtis data were further analyzed using the Adonis function in the Vegan package in R, and post-hoc testing was performed using the RVAideMemoire package using pairwise.perm.manova(), which conducts pairwise tests on matrix data using Adonis (Hervé, 2017; Oksanen et al., 2017).

reads<-as.data.frame(rowSums(tOTUList))
variables$reads<-reads

library("GUniFrac")
rareOTU<-Rarefy(tOTUList, depth = min(reads))  #Rarefaction
rffOTU<-as.data.frame(rareOTU$otu.tab.rff)  #2207
min(reads)
totalSeq<-rowSums(rffOTU)
trffOTU<-t(rffOTU)
raOTU <- rffOTU[ ,colSums(rffOTU)!=0]
traOTU<-t(raOTU)
paOTU<-raOTU
paOTU[paOTU>0]=1
row.names(raOTU)<-variables$uni

#Alpha diversity
richOTU<-rowSums(paOTU)
plot(richOTU~variables$rot_crus)
aovrich<-aov(log(richOTU)~variables$rot_crus)
tuk<-TukeyHSD(aovrich)
par(mfrow=c(1,2))
plot(richOTU~variables$env_cult)
aovrich<-aov(log(richOTU)~variables$env_cult)
tuk<-TukeyHSD(aovrich)

lmrich<-lm(richOTU~variables$rot_crus*variables$env_cult)

glm<-glm(richOTU~variables$env_cult*variables$substrate, family=quasipoisson(link="log"))

summary(glm)
env_cult<-variables$env_cult
rot_crus<-variables$rot_crus
glm<-glmer(richOTU~env_cult+rot_crus+(1|variables$Place), family=poisson)
glm2<-glmer(richOTU~rank+(1|variables$Place), family=poisson)
glm2<-glmer(richOTU~rank+(1|variables$Place), family=poisson)
Anova(glm)
library(multcomp)
summary(glht(glm, mcp(rot_crus="Tukey")))
summary(glht(glm2, mcp(rank="Tukey")))
Anova(glm2)
library(multcomp)
summary(glht(glm2, mcp(rank="Tukey")))

plot(richOTU~variables$ani_cult)
aovrich<-aov(log(richOTU)~variables$ani_cult) 
tuk<-TukeyHSD(aovrich)


shan<-diversity(raOTU, index="shannon")
plot(shan~variables$ani_cult)
aovrich<-aov(log(shan)~variables$ani_cult)
tuk<-TukeyHSD(aovrich)

evn<-shan/log(specnumber(raOTU))
plot(evn~variables$ani_cult)
aovrich<-aov(log(evn)~variables$ani_cult)
tuk<-TukeyHSD(aovrich)

k <- sample(nrow(raOTU), 1)
fish <- fisherfit(raOTU[k,])
plot(fish)

library("picante")

phy<-read.nexus("raxml_rootedotu.nexus")
comm<-raOTU
prunedphy <- prune.sample(comm, phy)
par(mfrow = c(2, 3))
pd.result <- pd(comm, phy, include.root=TRUE)
plot(pd.result$PD~variables$ani_cult)
par(mfrow=c(1,2))

rank<-variables$ani_cult
place<-variables$Place
glmR<-glm(pd.result$PD~rank, family=quasipoisson)
glmRich<-lmer(pd.result$PD~rank+(1|place)) #, family=gaussian) #Is there a statistical difference in richness between Microplastic and water samples?
summary(glmRich)
glmRichOUT<-capture.output(summary(glmRich)) #Write output in a way that is easily writable in a file
phydist <- cophenetic(phy)
#ses.mpd.result <- ses.mpd(comm, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
plot(richOTU~rank)
plot(ses.mpd.result$mpd.obs~rank)
glmRich<-glm(ses.mpd.result$mpd.obs~rank, family=quasipoisson) #Is there a statistical difference in richness between Microplastic and water samples?
summary(glmRich)
library(lmerTest)
lms<-lmer(ses.mpd.result$mpd.obs~rank+(1|variables$Place))
subs<-variables$substrate
glm<-lmer(ses.mpd.result$mpd.obs~env_cult+subs+(1|variables$Place)) #, family=poisson)
Anova(glm)

summary(lms)

library("car")
Alm<-Anova(lms)
Alm
Aglm<-Anova(glmRich)
Aglm
library(multcomp)
summary(glht(glmRich, mcp(rank="Tukey")))
summary(glht(glm2, mcp(rank="Tukey")))

glmRich<-glm(richOTU~rank, family=quasipoisson) #Is there a statistical difference in richness between Microplastic and water samples?
summary(glmRich)
summary(glht(glmRich, mcp(rank="Tukey")))

glmnb<-glm.nb(richOTU~rank)
library("performance")
check_model(glmnb)
check_model(glmRich)
