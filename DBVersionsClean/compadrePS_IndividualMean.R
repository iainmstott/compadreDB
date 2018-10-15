################################################
###CALCULATE PACE/SHAPE DATA FOR COMPADRE 3.0###
################################################
setwd("C:/Dropbox/Work/Data/Compadre")

###PACKAGES###
library(Rage)
library(ape)
library(popdemo)


###DATA###

#compadre version3.2.1, including phylogeny
load("phylo/compadre_v3.2.1-phy.Rdata")
#compadre version 3.2.1 plus Exeter extras
compadreplus<-read.csv("COMPADRE_v3.2.1-correctedIS.csv")

#check the two match
compadre$metadata$rowID<-apply(compadre$metadata[,c("SpeciesAuthor","Journal","YearPublication","LatNS","LonWE")],1,function(x){paste(as.character(x),collapse=" - ")})
compadreplus$rowID<-apply(compadreplus[,c("SpeciesAuthor", "Journal","YearPublication","LatNS","LonWE")],1,function(x){paste(as.character(x),collapse=" - ")})
all(compadre$metadata$rowID==compadreplus$rowID)

#migrate important information from compadreplus to compadre
compadre$metadata<-data.frame(compadre$metadata,IUCN=compadreplus$endangered,SuccessionGrowthType=compadreplus$growth_type_dave)

#work out which matrices contain NAs
compadre$metadata$hasna<-sapply(compadre$mat,function(x){any(is.na(x$matA))})
#work out which matrices are reducible
compadre$metadata$irr<-numeric(dim(compadre$metadata)[1]); compadre$metadata$irr[-which(compadre$metadata$hasna)]<-sapply(compadre$mat[-which(compadre$metadata$hasna)],function(x){is.matrix_irreducible(x$matA)}); compadre$metadata$irr<-as.logical(compadre$metadata$irr)
#work out which matrices have seeds problems
compadre$metadata$seedsp<-!is.na(compadreplus$Simons_Seeds_Corrected_pre_PPMS)&!is.na(compadreplus$Simons_Seeds_Corrected_post_PPMS)

#these are the matrices to use in the analysis...
#no NAs
#Irreducible
#No seeds problem
#In the phylogeny
#Annual, divided, unmanipulated, wild
use0<- !compadre$metadata$hasna & compadre$metadata$irr & !compadre$metadata$seedsp & compadre$metadata$InPhylo & compadre$metadata$MatrixSplit=="Divided" & compadre$metadata$AnnualPeriodicity=="1" & compadre$metadata$MatrixTreatment=="Unmanipulated" & compadre$metadata$MatrixCaptivity=="W"
#We only want to use individual matrices....
indM<-compadre$metadata$MatrixComposite=="Individual"
#...BUT we also want to use matrices that only have one entry per species
#and these could be recorded as many different things (Mean, Pooled, Individual)
numM<-table(compadre$metadata$SpeciesAccepted) #number of matrices per species
lonelySpp<-names(numM)[numM==1] #names of species with only one matrix
single<-compadre$metadata$SpeciesAccepted%in%lonelySpp #corresponding rows of metadata
#checks:
#1 are they really only single species?
length(which(single))==length(lonelySpp)
#2 what are the matrix types for the singles?
unique(compadre$metadata$MatrixComposite[single])
#no seasonals, good

#We want our matrices to be in use0, AND individual OR single
use1<- use0 & (indM | single)

#subset the data
compadre$metadata<-data.frame(compadre$metadata,use1)
compadre<-subsetDB(compadre,use1)
compadre$metadata<-compadre$metadata[,!names(compadre$metadata)=="use1"]



#Make life tables from these
maxage<-10000
lifetables<-lapply(compadre$mat,function(mat){makeLifeTable(matU=mat$matU,nSteps=maxage)})
#trim the life tables
zerocheck<-logical(length(lifetables))
for(i in 1:length(lifetables)){
    lifetables[[i]]$lx[lifetables[[i]]$lx<0.01]<-0           #any lx less than 0.01 becomes zero
    lifetables[[i]]$lx[maxage]<-0                            #if not then the last lx is always zero
}

#Run some checks on the life tables: 
#What are their dimensions? Should be greater than 3 and less than 10000
LTdim<-sapply(lifetables,function(LT){min(which(LT$lx==0))})
#check nothing comes back from the dead
for(i in 1:length(lifetables)){
    zerocheck[i]<-all(diff(which(lifetables[[i]]$lx==0))==1) #these are the ages (+1) which are zero: wanna check they're all consecutive (things don't arise from the dead)
}
#trim the life tables to do the rest
for(i in 1:length(lifetables)){
    lifetables[[i]]<-lifetables[[i]][1:LTdim[i],]
}
#are the life tables monotonically declining? (They should be)
LTmono<-sapply(lifetables,function(LT){all(diff(LT$lx)<0)})
#What's the max and min? (Should be 1 and 0 respectively),
LTmax<-sapply(lifetables,function(LT){max(LT$lx)})
LTmin<-sapply(lifetables,function(LT){min(LT$lx)})
#When is quasi convergence reached? Should be at age after LTdim
qsdc<-sapply(compadre$mat,function(mat){qsdConverge(mat$matU,nSteps=10000,conv=0.01)})
#Compile all this info together
LTok<- zerocheck & LTmono & !is.na(LTmax) & !is.na(LTmin) & LTmax==1 & LTmin==0 & !is.na(qsdc) & LTdim>3 & LTdim<maxage & LTdim<qsdc
#get rid of last class (where lx=0)
lifetables<-lapply(lifetables,function(x){x[1:(dim(x)[1]-1),]})
LTdim<-LTdim-1

#subset the data to those matrices with OK life tables
use2<-LTok
compadre$metadata<-data.frame(compadre$metadata,LTdim,qsdc,use2)
compadre<-subsetDB(compadre,use2)
compadre$metadata<-compadre$metadata[,!names(compadre$metadata)=="use2"]
lifetables<-lifetables[use2]



#calculate the pace measure: life expectancy
source("C:/Dropbox/Work/Software/lifehistory/lexp.R")
e0<-sapply(lifetables,function(lifetable){lexp(lifetable,x=0)})

#calculate the shape measure: Keyfitz's entropy
kent<-sapply(lifetables,function(lifetable){kentropy(lifetable$lx)})

#calculate lambdas
lambdamax<-sapply(compadre$mat,function(x){abs(Re(eigen(x$matA)$values[1]))})
#any extreme lambdas?
which(log(lambdamax)<(-1)|log(lambdamax)>1) #135, 137 and 152 are outside our preferred lambda range

#calculate transient bounds
reaclwr<-sapply(compadre$mat,function(x){firststepatt(x$matA)})
maxatt<-sapply(compadre$mat,function(x){maxatt(x$matA)})
inertlwr<-sapply(compadre$mat,function(x){inertia(x$matA,bound="lower")})
reacupr<-sapply(compadre$mat,function(x){reactivity(x$matA)})
maxamp<-sapply(compadre$mat,function(x){maxamp(x$matA)})
inertupr<-sapply(compadre$mat,function(x){inertia(x$matA,bound="upper")})

#Do principal components analysis on transients
pca<-princomp(log(data.frame(reaclwr,maxatt,inertlwr,reacupr,maxamp,inertupr)))
par(mfrow=c(1,2)); plot(pca); biplot(pca) #151 is an outlier

#take out the outliers
use3<-logical(length(lifetables)); use3[-c(135,137,151,152)]<-TRUE
compadre$metadata<-data.frame(compadre$metadata,use3)
compadre<-subsetDB(compadre,use3)
compadre$metadata<-compadre$metadata[,!names(compadre$metadata)=="use3"]
lifetables<-lifetables[use3]
e0<-e0[use3]
kent<-kent[use3]
lambdamax<-lambdamax[use3]
reaclwr<-reaclwr[use3]
maxatt<-maxatt[use3]
inertlwr<-inertlwr[use3]
reacupr<-reacupr[use3]
maxamp<-maxamp[use3]
inertupr<-inertupr[use3]

#Check lambdamax again
which(log(lambdamax)<(-1)|log(lambdamax)>1) #no outliers
#Do principal components analysis on transients again
pca<-princomp(log(data.frame(reaclwr,maxatt,inertlwr,reacupr,maxamp,inertupr)))
par(mfrow=c(1,2)); plot(pca); biplot(pca) #no outliers

#Add to data (only PC1 and PC2 needed)
compadre$metadata<-data.frame(compadre$metadata,e0,kent,lambdamax,reaclwr,maxatt,inertlwr,reacupr,maxamp,inertupr,PC1=pca$score[,1],PC2=pca$score[,2])



allvars<-c("SpeciesAuthor","SpeciesAccepted","Authority","TaxonomicStatus","TPLVersion",
          "InfraspecificAccepted","SpeciesEpithetAccepted","GenusAccepted","Genus",
          "Family","Order","Class","DicotMonoc","AngioGymno","Phylum","Kingdom",
          "MatrixDimension","GrowthType","InPhylo","IUCN","SuccessionGrowthType")
PS_IndMean<-compadre$metadata[!duplicated(compadre$metadata$SpeciesAccepted),allvars]
modelvars<-t(sapply(PS_IndMean$SpeciesAccepted,function(x){colMeans(compadre$metadata[compadre$metadata$SpeciesAccepted==x,c("MatrixDimension","e0","kent","lambdamax","reaclwr","maxatt","inertlwr","reacupr","maxamp","inertupr")])}))
PS_IndMean<-data.frame(PS_IndMean,modelvars)

#Do principal components analysis on mean transients
pcaMean<-princomp(log(PS_IndMean[,c("reaclwr","maxatt","inertlwr","reacupr","maxamp","inertupr")]))
par(mfrow=c(1,2)); plot(pcaMean); biplot(pcaMean) #no outliers

#Add to mean data (only PC1 and PC2 needed)
PS_IndMean<-data.frame(PS_IndMean,PC1=pcaMean$score[,1],PC2=pcaMean$score[,2])



#create new objects
compadre$metadata<-compadre$metadata[,!names(compadre$metadata)%in%c("rowID")]
phy<-drop.tip(phy,phy$tip.label[!phy$tip.label%in%compadre$metadata$SpeciesAccepted])
all(phy$tip.label%in%compadre$metadata$SpeciesAccepted)
all(compadre$metadata$SpeciesAccepted%in%phy$tip.label)
obj<-ls()
rm(list=obj[!obj=="compadre"&!obj=="phy"&!obj=="lifetables"&!obj=="pca"&!obj=="PS_IndMean"&!obj=="pcaMean"]); rm(obj)
save.image("PaceShape/compadrePS_Individual.Rdata")


