setwd("C:/Users/Ben/Dropbox/Malaria/Data_analysis/Parity")
rm(list=ls())
library(rstan)
library(shinystan)
options(mc.cores=parallel::detectCores())
aDF <- read.csv('aggregatedReordered_v1.csv')


## Removing series where there is insecticide, or removing the nulliparous observations if
## Np < X * P1, or removing series with insufficient observations
fPrepareSpecies <- function(aDF,aSpeciesNumber,aThreshold,aThresholdSize){
  aDF_gambiae <- aDF[aDF$SpeciesID==aSpeciesNumber,]
  aDF_gambiae <- aDF_gambiae[aDF_gambiae$Control==1,]
  aDF_gambiae <- aDF_gambiae[aDF_gambiae$Total>=aThresholdSize,]
  lRemoveNP <- aDF_gambiae$Nulliparous<aThreshold*aDF_gambiae$X1P
  obsDF <- aDF_gambiae[,c(11:24,27)]
  lThreshold <- ifelse(obsDF$Threshold>0,obsDF$Threshold,-1)
  lStartingT <- ifelse(lRemoveNP,1,0)
  
  ## If there are no studies left then return -1
  if(nrow(aDF_gambiae)==0){
    return(list(K=-1))
  }
  
  ## Check whether thresholding has been recorded corrected
  K <- nrow(obsDF)
  for(i in 1:K){
    aTempDF <- obsDF[i,]
    if(lThreshold[i]>0){
      if(lThreshold[i] < 13){
        if(sum(aTempDF[(lThreshold[i]+1):13])>0){
          print(paste0('an error occurred in ...',i))
        }
      }
      
    }
  }
  
  mObs <- matrix(nrow = K,ncol = 14)
  for(i in 1:K){
    mObs[i,1:13] <- unlist(obsDF[i,1:13])
    mObs[i,14] <- 0
    if(lThreshold[i]>0){
      mObs[i,(lThreshold[i]+1)] <- obsDF[i,14]
    }
  }
  
  if(sum(mObs) != sum(obsDF[,c(1:14)]))
    print('an error with sums...')
  
  ## stack capture observations
  lY <- vector()
  lS <- vector(length = K)
  for(i in 1:K){
    if(!lRemoveNP[i]){
      if(lThreshold[i]>0){
        lY <- c(lY,mObs[i,1:(lThreshold[i]+1)])
        lS[i] <- length(mObs[i,1:(lThreshold[i]+1)])
      } else{
        lY <- c(lY,mObs[i,])
        lS[i] <- length(mObs[i,])
      }
    }else{
      if(lThreshold[i]>0){
        lY <- c(lY,mObs[i,2:(lThreshold[i]+1)])
        lS[i] <- length(mObs[i,2:(lThreshold[i]+1)])
      } else{
        aLen <- length(mObs[i,])
        lY <- c(lY,mObs[i,2:aLen])
        lS[i] <- length(mObs[i,2:aLen])
      }
    }
  }
  
  if(length(lY) != sum(lS))
    print('an error with lengths')
  N <- length(lY)
  Pos <- c(1,(cumsum(lS)+1))
  Pos <- Pos[1:(length(Pos)-1)]
  ID=as.character(aDF_gambiae$ID....concat)
  
  ## Check sums
  for(i in 1:K){
    aSumA <- sum(lY[Pos[i]:(Pos[i]+lS[i]-1)])
    aSumB <- ifelse(!lRemoveNP[i],
                    aDF_gambiae$Total[aDF_gambiae$ID....concat==ID[i]],
                    aDF_gambiae$Total[aDF_gambiae$ID....concat==ID[i]]-aDF_gambiae$Nulliparous[aDF_gambiae$ID....concat==ID[i]])
    if(aSumA!=aSumB)
      print(paste0('an error in i = ',i))
  }
  
  return(list(N=N,Y=lY,Pos=as.array(Pos),K=K,S=as.array(lS),
              threshold=as.array(lThreshold),removeNP=lRemoveNP,ID=ID,
              startingT=lStartingT))
}


aDF1 <- fPrepareSpecies(aDF,1,0.9,100)

dModel <- stan_model('dissectionNB1_exponential.stan')
initFnGompertz <- function(){
  list(alpha=runif(K,min = 0,max=0.6),
       beta=runif(K,min = 1,max=1.6))
}
aTempDF <- fPrepareSpecies(aDF,1,0.9,100)
K <- aTempDF$K
fit <- sampling(dModel,data=fPrepareSpecies(aDF,8,0.9,100),iter=200,chains=4)
print(fit)
lAlpha <- extract(fit,'alpha_average')[[1]]
lBeta <- extract(fit,'beta_average')[[1]]
lAlphaBeta <- data.frame(alpha=lAlpha,beta=lBeta)
curve(exp(-0.15/0.33 * exp(0.33*x)),0,10)

aIntegralUpper <- 1000
fCalculateLifetime_gompertz <- function(overallAlpha,overallBeta){
  aInt <- integrate(function(t) exp(-(overallAlpha/overallBeta)*(exp(overallBeta*t)-1)),0,aIntegralUpper)
  return(aInt[[1]])
}

lLifetimes <- apply(lAlphaBeta,1,function(x) 
  do.call(fCalculateLifetime_gompertz,unname(as.list(x))))
quantile(lLifetimes)

fCalculateLifetime_gompertzMakeham <- function(overallAlpha,overallBeta,overallcMakeham){
  aInt <- integrate(function(t) exp(-(overallAlpha/overallBeta)*(exp(overallBeta*t)-1)-overallcMakeham*t),0,aIntegralUpper)
  return(aInt[[1]])
}

fCalculateLifetime_logistic <- function(overallAlpha,overallBeta,overallSS){
  lInt=try(integrate(function(t) (((overallBeta-overallAlpha*overallSS+
                                      overallAlpha*overallSS*exp(overallBeta*t))/overallBeta)^(-1/overallSS)),0,aIntegralUpper),
           TRUE)
  if(is.atomic(lInt)){lInt=0; return(lInt)}
  return(lInt$value)
}

fCalculateLifetime_logisticMakeham <- function(overallAlpha,overallBeta,overallSS,overallcMakeham){
  lInt=try(integrate(function(t) exp(-overallcMakeham*t)*(((overallBeta-overallAlpha*overallSS+
                                                              overallAlpha*overallSS*exp(overallBeta*t))/overallBeta)^(-1/overallSS)),0,aIntegralUpper),
           TRUE)
  if(is.atomic(lInt)){lInt=0; return(lInt)}
  return(lInt$value)
}

fCalculateLifetime_weibull <- function(overallsambda,overallBeta){
  return(((overallsambda/overallBeta)^(-1/overallBeta)) * gamma(1+(1/overallBeta)))
}


lAlpha <- extract(fit,'alpha_average')[[1]]
lBeta <- extract(fit,'beta_average')[[1]]
lSS <- extract(fit,'ss_average')[[1]]
lcMakeham <- extract(fit,'cMakeham_average')[[1]]
lAlphaBeta <- data.frame(alpha=lAlpha,beta=lBeta,ss=lSS,cMakeham=lcMakeham)

lLifetimes <- apply(lAlphaBeta,1,function(x) 
  do.call(fCalculateLifetime_logisticMakeham,unname(as.list(x))))
quantile(lLifetimes)





## Np < X * P1, or removing series with insufficient observations
fPrepareSpeciesHoldOut <- function(aDF,aSpeciesNumber,aThreshold,aThresholdSize,
                                   aIncludeFold,aHoldOutFold){
  aDF_gambiae <- aDF[aDF$SpeciesID==aSpeciesNumber,]
  aDF_gambiae <- aDF_gambiae[aDF_gambiae$Control==1,]
  aDF_gambiae <- aDF_gambiae[aDF_gambiae$Total>=aThresholdSize,]
  lRemoveNP <- aDF_gambiae$Nulliparous<aThreshold*aDF_gambiae$X1P
  obsDF <- aDF_gambiae[,c(11:24,27)]
  lThreshold <- ifelse(obsDF$Threshold>0,obsDF$Threshold,-1)
  lStartingT <- ifelse(lRemoveNP,1,0)
  
  ## If there are no studies left then return -1
  if(nrow(aDF_gambiae)==0){
    return(list(K=-1))
  }
  
  ## Check whether thresholding has been recorded corrected
  K <- nrow(obsDF)
  for(i in 1:K){
    aTempDF <- obsDF[i,]
    if(lThreshold[i]>0){
      if(lThreshold[i] < 13){
        if(sum(aTempDF[(lThreshold[i]+1):13])>0){
          print(paste0('an error occurred in ...',i))
        }
      }
      
    }
  }
  
  mObs <- matrix(nrow = K,ncol = 14)
  for(i in 1:K){
    mObs[i,1:13] <- unlist(obsDF[i,1:13])
    mObs[i,14] <- 0
    if(lThreshold[i]>0){
      mObs[i,(lThreshold[i]+1)] <- obsDF[i,14]
    }
  }
  
  if(sum(mObs) != sum(obsDF[,c(1:14)]))
    print('an error with sums...')
  
  ## stack capture observations
  lY <- vector()
  lS <- vector(length = K)
  for(i in 1:K){
    if(!lRemoveNP[i]){
      if(lThreshold[i]>0){
        lY <- c(lY,mObs[i,1:(lThreshold[i]+1)])
        lS[i] <- length(mObs[i,1:(lThreshold[i]+1)])
      } else{
        lY <- c(lY,mObs[i,])
        lS[i] <- length(mObs[i,])
      }
    }else{
      if(lThreshold[i]>0){
        lY <- c(lY,mObs[i,2:(lThreshold[i]+1)])
        lS[i] <- length(mObs[i,2:(lThreshold[i]+1)])
      } else{
        aLen <- length(mObs[i,])
        lY <- c(lY,mObs[i,2:aLen])
        lS[i] <- length(mObs[i,2:aLen])
      }
    }
  }
  
  if(length(lY) != sum(lS))
    print('an error with lengths')
  N <- length(lY)
  Pos <- c(1,(cumsum(lS)+1))
  Pos <- Pos[1:(length(Pos)-1)]
  ID=as.character(aDF_gambiae$ID....concat)
  
  ## Check sums
  for(i in 1:K){
    aSumA <- sum(lY[Pos[i]:(Pos[i]+lS[i]-1)])
    aSumB <- ifelse(!lRemoveNP[i],
                    aDF_gambiae$Total[aDF_gambiae$ID....concat==ID[i]],
                    aDF_gambiae$Total[aDF_gambiae$ID....concat==ID[i]]-aDF_gambiae$Nulliparous[aDF_gambiae$ID....concat==ID[i]])
    if(aSumA!=aSumB)
      print(paste0('an error in i = ',i))
  }
  
  K2 <- length(aHoldOutFold)
  N2 <- sum(lS[aHoldOutFold])-K2;
  
  return(list(N=N,Y=lY,Pos=as.array(Pos),K1=(K-K2),K=K,S=as.array(lS),
              threshold=as.array(lThreshold),removeNP=lRemoveNP,ID=ID,
              startingT=lStartingT,includeSeries=aIncludeFold,holdOutSeries=aHoldOutFold,
              K2=K2,N2=N2))
}

library(caret)
lK <- vector(length=length(unique(aDF$Species)))
for(i in 1:length(unique(aDF$Species))){
  print(i)
  aTempDF <- fPrepareSpecies(aDF,i,0.9,100)
  lK[i] <- aTempDF$K
}

lSufficient <- which(lK>1)
lFolds <- vector(length = length(lSufficient),mode='list')
for(i in seq_along(lSufficient)){
  lFolds[[i]] <- createFolds(1:lK[lSufficient[i]],2)
}
save(file = 'lFolds.RDa',lFolds)

load('lFolds.RDa')

options(mc.cores = 4)

bModel <- stan_model('dissectionNB1_gompertz_kFold.stan')
aModel <- stan_model('dissectionNB1_exponential_kFold.stan')
cModel <- stan_model('dissectionNB1_logisticMakeham_kFold.stan')
aSpecies <- 1
fit_gompertz <- sampling(cModel,data=fPrepareSpeciesHoldOut(aDF,aSpecies,0.9,100,lFolds[[aSpecies]]$Fold1,lFolds[[aSpecies]]$Fold2),
                    iter=200,chains=4)
fit_exp <- sampling(aModel,data=fPrepareSpeciesHoldOut(aDF,aSpecies,0.9,100,lFolds[[aSpecies]]$Fold1,lFolds[[aSpecies]]$Fold2),
                iter=200,chains=4)
library(loo)
library(matrixStats)
lLog_exp <- extract_log_lik(fit1,'lLogLikelihood')
lLog_gompertz <- extract_log_lik(fit_gompertz,'lLogLikelihood')



aLen <- dim(lLog_exp)[2]
bLen <- dim(lLog_exp)[1]
lSumExp <- sapply(seq(1,aLen,1),function(i) logSumExp(lLog_gompertz[,i])) - log(bLen)
sum(lSumExp)

fGetlpd <- function(lMLogProb){
  
  aLP <- lMLogProb
  numFolds <- length(aLP)
  lLPD <- vector()
  for (j in 1:numFolds){
    aTest <- aLP[[j]]
    aLen <- dim(aTest)[2]
    bLen <- dim(aTest)[1]
    print(aLen)
    lSumExp <- sapply(seq(1,aLen,1),function(i) logSumExp(aTest[,i])) - log(bLen)
    lLPD <- c(lLPD,lSumExp)
  }
  return(lLPD)
}


fPrepareAll <- function(aDF,aThreshold,aThresholdSize){
  aDF_gambiae <- aDF[aDF$Control==1,]
  aDF_gambiae <- aDF_gambiae[aDF_gambiae$Total>=aThresholdSize,]
  aTable <- as.data.frame(table(aDF_gambiae$SpeciesID))
  lSpeciesRemove <- aTable$Var1[aTable$Freq<2]
  aDF_gambiae <- aDF_gambiae[!aDF_gambiae$SpeciesID%in%lSpeciesRemove,]
  lRemoveNP <- aDF_gambiae$Nulliparous<aThreshold*aDF_gambiae$X1P
  obsDF <- aDF_gambiae[,c(11:24,27)]
  lThreshold <- ifelse(obsDF$Threshold>0,obsDF$Threshold,-1)
  lStartingT <- ifelse(lRemoveNP,1,0)
  lSpeciesID <- aDF_gambiae$SpeciesID
  aLenSpecies <- length(unique(lSpeciesID))
  lGenus <- aDF_gambiae$Genus
  aLenGenus <- length(unique(lGenus))
  
  cLookupGenus <- data.frame(genus_name=unique(lGenus),genus=1:aLenGenus)
  lGenus_new <- cLookupGenus$genus[match(lGenus,cLookupGenus$genus_name)]
  
  bLookup <- data.frame(newSpecies=1:aLenSpecies,old=unique(lSpeciesID))
  lSpeciesID_new <- bLookup$newSpecies[match(lSpeciesID,bLookup$old)]
  
  ## If there are no studies left then return -1
  if(nrow(aDF_gambiae)==0){
    return(list(K=-1))
  }
  
  ## Check whether thresholding has been recorded corrected
  K <- nrow(obsDF)
  for(i in 1:K){
    aTempDF <- obsDF[i,]
    if(lThreshold[i]>0){
      if(lThreshold[i] < 13){
        if(sum(aTempDF[(lThreshold[i]+1):13])>0){
          print(paste0('an error occurred in ...',i))
        }
      }
      
    }
  }
  
  mObs <- matrix(nrow = K,ncol = 14)
  for(i in 1:K){
    mObs[i,1:13] <- unlist(obsDF[i,1:13])
    mObs[i,14] <- 0
    if(lThreshold[i]>0){
      mObs[i,(lThreshold[i]+1)] <- obsDF[i,14]
    }
  }
  
  if(sum(mObs) != sum(obsDF[,c(1:14)]))
    print('an error with sums...')
  
  ## stack capture observations
  lY <- vector()
  lS <- vector(length = K)
  for(i in 1:K){
    if(!lRemoveNP[i]){
      if(lThreshold[i]>0){
        lY <- c(lY,mObs[i,1:(lThreshold[i]+1)])
        lS[i] <- length(mObs[i,1:(lThreshold[i]+1)])
      } else{
        lY <- c(lY,mObs[i,])
        lS[i] <- length(mObs[i,])
      }
    }else{
      if(lThreshold[i]>0){
        lY <- c(lY,mObs[i,2:(lThreshold[i]+1)])
        lS[i] <- length(mObs[i,2:(lThreshold[i]+1)])
      } else{
        aLen <- length(mObs[i,])
        lY <- c(lY,mObs[i,2:aLen])
        lS[i] <- length(mObs[i,2:aLen])
      }
    }
  }
  
  if(length(lY) != sum(lS))
    print('an error with lengths')
  N <- length(lY)
  Pos <- c(1,(cumsum(lS)+1))
  Pos <- Pos[1:(length(Pos)-1)]
  ID=as.character(aDF_gambiae$ID....concat)
  
  ## Check sums
  for(i in 1:K){
    aSumA <- sum(lY[Pos[i]:(Pos[i]+lS[i]-1)])
    aSumB <- ifelse(!lRemoveNP[i],
                    aDF_gambiae$Total[aDF_gambiae$ID....concat==ID[i]],
                    aDF_gambiae$Total[aDF_gambiae$ID....concat==ID[i]]-aDF_gambiae$Nulliparous[aDF_gambiae$ID....concat==ID[i]])
    if(aSumA!=aSumB)
      print(paste0('an error in i = ',i))
  }
  
  return(list(N=N,Y=lY,Pos=as.array(Pos),K=K,S=as.array(lS),
              threshold=as.array(lThreshold),removeNP=lRemoveNP,ID=ID,
              startingT=lStartingT,species=lSpeciesID_new,species_old=lSpeciesID,
              genus=lGenus_new,genus_old=lGenus,nGenera=aLenGenus,
              nSpecies=aLenSpecies))
}

lTest <- fPrepareAll(aDF,0.9,100)

# look at years
ids_1 <- lTest$ID
testy <- aDF %>% filter(ID....concat%in%ids_1)
testy$End.date <- as.character(testy$End.date)
testy <- testy %>% 
  mutate(year=substr(Start.date, 7, 11))
mean(testy$year < 2000, na.rm = T)

cModel <- stan_model('dissectionNB1_exponential_species_single.stan')
options(mc.cores=16)
fit1 <- sampling(cModel,data=lTest,iter=200,chains=16)
lAlpha_indiv <- extract(fit1,'alpha')[[1]]
lLifetimes_indiv <- 1/lAlpha_indiv
write.csv(file='lifetimes_individual_indiv.csv',lLifetimes_indiv)
lLifetimes_indiv[1]
lTemp <- read.csv('lifetimes_individual_hierarchical.csv')
lTemp[1,2]

print(fit1)
lAlpha <- extract(fit1,'alpha_average')[[1]]
lLifetimes <- 1/lAlpha
write.csv(file = 'lifetimes_exponential_species.csv',lLifetimes)
aLookup <- unique(data.frame(speciesID=aDF$SpeciesID,speciesName=aDF$Species,genus=aDF$Genus))
lTest$SpeciesName <- aLookup$speciesName[match(lTest$species_old,aLookup$speciesID)]
lTest$GenusName <- aLookup$genus[match(lTest$species_old,aLookup$speciesID)]
aTable <- as.data.frame(table(lTest$SpeciesName))
aTable <- aTable[aTable$Freq>0,]
bLookup <- unique(data.frame(speciesName=lTest$SpeciesName,genusName=lTest$GenusName))
aTable$Genus <- bLookup$genusName[match(aTable$Var1,bLookup$speciesName)]
bLookup$number <- aTable$Freq[match(bLookup$speciesName,aTable$Var1)]
write.csv(file = 'genusSpeciesNumbers.csv',bLookup)

aModel <- stan_model('dissectionNB1_exponential_genus.stan')
fit2 <- sampling(aModel,data=lTest,iter=200,chains=16)
print(fit2)
lLifetimes_genus <- 1/extract(fit2,'alpha_average')[[1]]
write.csv(file = 'lifetimes_exponential_genus.csv',lLifetimes_genus)
bTable <- as.data.frame(table(lTest$genus))
lTemp <- bTable$Freq
write.csv(file = 'lNumberPerGenus.csv',lTemp)

bModel <- stan_model('dissectionNB1_exponential_overall.stan')
fit3 <- sampling(bModel,data=lTest,iter=200,chains=16)
print(fit3)
lLifetimes_overall <- 1/extract(fit3,'alpha_average')[[1]]
write.csv(file = 'lifetimes_exponential_overall.csv',lLifetimes_overall)

eModel <- stan_model('dissectionNB1_exponential_individuals.stan')

options(mc.cores=16)
fit <- sampling(eModel,data=lTest,iter=200,chains=16)
print(fit)
lLifetimes_indiv <- 1/extract(fit,'alpha')[[1]]


mTest <- lLifetimes_indiv
lQuantiles <- apply(mTest,2,function(x) c(low=quantile(x,0.25,na.rm=T),med=median(x,na.rm=T),
                                          high=quantile(x,0.75,na.rm=T)))
lQuantiles <- t(as.data.frame(lQuantiles))
colnames(lQuantiles) <- c("high","med","low")
lQuantiles <- as.data.frame(lQuantiles)

lQuantiles$id <- seq(1,nrow(lQuantiles),1)
limits <- aes(ymax=high,ymin=low,colour=Genus)
lSpecies <- aLookup$speciesName[match(lTest$species_old,aLookup$speciesID)]
lGenus <- aLookup$genus[match(lTest$species_old,aLookup$speciesID)]
lQuantiles$Species <- as.factor(lSpecies)
lQuantiles$Genus <- lGenus


## Order estimates according to the median posterior estimate of each species
lMedian <- vector(length=nrow(aTable))
for(i in 1:nrow(aTable)){
  aTempDF <- lQuantiles[lQuantiles$Species==aLookup$speciesName[i],]
  lMedian[i] <- median(aTempDF$med)
}
lOrder <- as.data.frame(order(lMedian))
lOrder$new <- seq(1,nrow(aTable),1)
colnames(lOrder)[1] <- "old"

lQuantiles1 <- lQuantiles[0,]
for(i in 1:nrow(aTable)){
  aSpecies <- lOrder$old[i]
  aTempDF <- lQuantiles[lQuantiles$Species==aLookup$speciesName[aSpecies],]
  lQuantiles1 <- rbind.data.frame(lQuantiles1,aTempDF)
}

## Sort within species
lQuantiles2 <- lQuantiles[0,]
for(i in 1:nrow(aTable)){
  aSpecies <- lOrder$old[i]
  aTempDF <- lQuantiles[lQuantiles$Species==aLookup$speciesName[aSpecies],]
  aTempDF <- aTempDF[order(aTempDF$med),]
  lQuantiles2 <- rbind.data.frame(lQuantiles2,aTempDF)
}
lQuantiles1 <- lQuantiles2

lQuantiles1$x <- 0
for(i in 2:nrow(lQuantiles)){
  if(lQuantiles1$Species[i]==lQuantiles1$Species[i-1]){
    lQuantiles1$x[i] <- lQuantiles1$x[i-1] + 3
  } else{
    lQuantiles1$x[i] <- lQuantiles1$x[i-1] + 6
  }
}

lQuantiles1$number <- seq(1,nrow(aTable),1)
lLabelNumber <- vector(length=)
lLines <- vector(length=39)

for(i in 1:39){
  aSpecies <- lOrder$old[i]
  aTempDF <- lQuantiles1[lQuantiles1$Species==aSpecies,]
  aMax <- max(aTempDF$x)
  aMin <- min(aTempDF$x)
  lLabelNumber[i] <- (aMin+aMax)/2
  lLines[i] <- aMax + 3
}
lLines <- c(-3,lLines)

speciesName <- read.csv('speciesNameLookup.csv')$SpeciesName
speciesName <- sapply(speciesName,as.character)
speciesName <- speciesName[lOrder$old]
# lLabelNumber <- sapply(lLabelNumber,as.character)

lSpeciesName <- vector()
for(i in 1:39){
  aSpecies <- lOrder$old[i]
  aTempDF <- lQuantiles1[lQuantiles1$Species==aSpecies,]
  lSpeciesName <- c(lSpeciesName,rep(speciesName[i],dim(aTempDF)[1]))
}
lQuantiles1$SpeciesName <- lSpeciesName

Palette1 <-c(rgb(31/255,119/255,180/255),rgb(255/255,127/255,14/255),rgb(44/255,160/255,44/255))
g1<-ggplot(lQuantiles1,aes(x=x,y=med,colour=as.factor(Genus)))+geom_point(aes(colour=as.factor(Genus)),size=2)+
  scale_x_continuous(breaks=lLabelNumber,labels=speciesName) +
  geom_errorbar(limits)  + coord_cartesian(ylim=c(0,25))  + xlab('species') +
  ylab('mean lifespan estimate, days') + theme_bw(base_size = 24)+
  theme(panel.grid.major = element_line(size = .5, color = "grey"),
        axis.line = element_line(size=.7, color = "black"),
        axis.text.y = element_text(face="italic"),
        legend.text = element_text(face="italic",size=24),
        legend.title = element_text(size=24),
        legend.position = c(.85,.2),
        #increase the font size
        text = element_text(size=24)) + coord_flip() + scale_colour_manual(values=Palette1,name="Genus") +
  guides(colour = guide_legend(override.aes = list(size=10,linetype=0)))

aLen <- dim(lQuantiles1)[1]
lQuantiles1 <- lQuantiles1[-((aLen-1):aLen),]

g1<-ggplot(lQuantiles1,aes(x=x,y=med,colour=as.factor(Genus)))+geom_point(aes(colour=as.factor(Genus)),size=2)+
  scale_x_continuous(breaks=lLabelNumber[-c(39)],labels=speciesName[-c(39)]) +
  geom_errorbar(limits)  + coord_cartesian(ylim=c(0,25))  + xlab('species') +
  ylab('mean lifespan estimate, days') + theme_bw(base_size = 24)+
  theme(panel.grid.major = element_line(size = .5, color = "grey"),
        axis.line = element_line(size=.7, color = "black"),
        axis.text.y = element_text(face="italic"),
        legend.text = element_text(face="italic",size=24),
        legend.title = element_text(size=24),
        legend.position = c(.85,.2),
        #increase the font size
        text = element_text(size=24)) + coord_flip() + scale_colour_manual(values=Palette1,name="Genus") +
  guides(colour = guide_legend(override.aes = list(size=10,linetype=0)))


# Trying to reformat the graphs so that lines occur not in breaks
bLen <- length(lLines)
g1<-ggplot(lQuantiles1,aes(x=x,y=med,colour=as.factor(Genus)))+geom_point(aes(colour=as.factor(Genus)),size=2)+
  scale_x_continuous(breaks=lLabelNumber,labels=speciesName) +
  geom_errorbar(limits)   + xlab('species') +
  ylab('mean lifespan estimate, days') + theme_bw(base_size = 24)+
  theme(panel.grid.major = element_line(size = 0.0, color = "white"),
        panel.grid.minor= element_line(size = 0.0, color = "white"),
        axis.line = element_line(size=.7, color = "black"),
        axis.text.y = element_text(face="italic"),
        legend.text = element_text(face="italic",size=24),
        legend.title = element_text(size=24),
        legend.position = c(.85,.2),
        #increase the font size
        text = element_text(size=24)) + coord_flip() + scale_colour_manual(values=Palette1,name="Genus") +
  guides(colour = guide_legend(override.aes = list(size=10,linetype=0))) + geom_vline(xintercept = lLines,color="grey")

ggsave('individualEstimates_allSpecies.png',width= 24, height = 24,g1,dpi=400)
ggsave('individualEstimates_allSpecies.pdf',width= 24, height = 24,g1)

aLen <- dim(lQuantiles1)[1]
lQuantiles1 <- lQuantiles1[-((aLen-1):aLen),]


g1<-ggplot(lQuantiles1,aes(x=x,y=med,colour=as.factor(Genus)))+geom_point(aes(colour=as.factor(Genus)),size=2)+
  scale_x_continuous(breaks=lLabelNumber[-c(39)],labels=speciesName[-c(39)]) +
  geom_errorbar(limits)   + xlab('species') +
  ylab('mean lifespan estimate, days') + theme_bw(base_size = 24)+
  theme(panel.grid.major = element_line(size = 0.0, color = "white"),
        panel.grid.minor= element_line(size = 0.0, color = "white"),
        axis.line = element_line(size=.7, color = "black"),
        axis.text.y = element_text(face="italic"),
        legend.text = element_text(face="italic",size=24),
        legend.title = element_text(size=24),
        legend.position = c(.85,.2),
        #increase the font size
        text = element_text(size=24)) + coord_flip() + scale_colour_manual(values=Palette1,name="Genus") +
  guides(colour = guide_legend(override.aes = list(size=10,linetype=0))) + geom_vline(xintercept = lLines[-c(bLen)],color="grey")

ggsave('individualEstimates_allSpecies_withoutBalabacensis.png',width= 24, height = 24,g1,dpi=400)
ggsave('individualEstimates_allSpecies_withoutBalabacensis.pdf',width= 24, height = 24,g1,dpi=400)

# Log x axis 
bLen <- length(lLines)
g1<-ggplot(lQuantiles1,aes(x=x,y=med,colour=as.factor(Genus)))+geom_point(aes(colour=as.factor(Genus)),size=2)+
  scale_x_continuous(breaks=lLabelNumber,labels=speciesName) + scale_y_log10()+
  geom_errorbar(limits)   + xlab('species') +
  ylab('mean lifespan estimate, days') + theme_bw(base_size = 24)+
  theme(panel.grid.major = element_line(size = 0.0, color = "white"),
        panel.grid.minor= element_line(size = 0.0, color = "white"),
        axis.line = element_line(size=.7, color = "black"),
        axis.text.y = element_text(face="italic"),
        legend.text = element_text(face="italic",size=24),
        legend.title = element_text(size=24),
        legend.position = c(.85,.2),
        #increase the font size
        text = element_text(size=24)) + coord_flip() + scale_colour_manual(values=Palette1,name="Genus") +
  guides(colour = guide_legend(override.aes = list(size=10,linetype=0))) + geom_vline(xintercept = lLines,color="grey")

ggsave('individualEstimates_allSpecies_logScale.pdf',width= 24, height = 24,g1,dpi=400)



