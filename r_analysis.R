#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#read in metadata and asv table
set.seed(515)
asv_table <- read.delim("~/Documents/GitHub/JuliaCircadian/julia_asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/Documents/GitHub/JuliaCircadian/julia_synthcom_daynight_map.txt", header=T)


#look at sequencing depth
min(colSums(asv_table))
#2635

max(colSums(asv_table))
#22186

#rarefy data 
nut_rare<-rrarefy(t(asv_table), sample=2635)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(nut_rare  ~ 1, distance='bray')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#36.4
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#14.8

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color=Type))+
  geom_point()+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 36.4%")+
  ylab("PC2- 14.8%")



###Calculate alpha diversity


#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(asv_table), sample=2635)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')

larv.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(asv_table), sample=2635), index = 'shannon'))
names(larv.alpha2)<-"Shannon"
larv.alph<-cbind(larv.alph, larv.alpha2)


#plot richness and shannon
ggplot(larv.alph, aes(Type, Shannon, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")

ggplot(larv.alph, aes(Type, `specnumber(rrarefy(t(asv_table), sample = 2635))`, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("sOTUs Observed")


