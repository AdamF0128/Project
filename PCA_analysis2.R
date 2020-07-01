#LOADING THE DATA
rm(list=ls())
#dev.off(dev.list()["RStudioGD"]) #Clear all plots

#These lines load the required library files into the workspace
#If there is an error here then use install.packages("package name",dependencies=TRUE)
#at the command line to perform an install
library(R.matlab)
library(hyperSpec)
library(ChemometricsWithR)
library(baseline)
library(ggplot2)
library(pls)

#Loading wavenumber data
wavenumber=readMat("wavenumber.mat")
wavenumber=wavenumber$wavenumber
###################LNCaP#######################
#LNCaP
data=readMat("LNCaP_0Gy.mat") #This loads the matlab file into the workspace as a 'list' object
LNCaP.0Gy=data$data #Extracting the LNCaP spectral data matrix

data=readMat("LNCaP_IF.mat") #This loads the matlab file into the workspace as a 'list' object
LNCaP.IF=data$data 

#SG Filtering
numrow=nrow(LNCaP.0Gy)
numcol=ncol(LNCaP.0Gy)
vec1=NULL
vec1[seq(1,numrow)]="LNCaP 0Gy"

filtLNCaP.0Gy=as.matrix(LNCaP.0Gy,nrow=numrow,ncol=numcol)
library(signal)

for (i in 1:numrow){
  filtLNCaP.0Gy[i,]=sgolayfilt(LNCaP.0Gy[i,],p=5,n=13,m=0)
}

print(LNCaP.0Gy)
#Baseline Correction
LNCaP.0Gy <- baseline(filtLNCaP.0Gy, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
LNCaP.0Gy <- LNCaP.0Gy@corrected

#SG Filtering
numrow1=nrow(LNCaP.IF)
numcol1=ncol(LNCaP.IF)
vec2=NULL
vec2[seq(1,numrow1)]="LNCaP 2Gy"

filtLNCaP.IF=as.matrix(LNCaP.IF,nrow=numrow1,ncol=numcol1)

for (i in 1:numrow1){
  filtLNCaP.IF[i,]=sgolayfilt(LNCaP.IF[i,],p=5,n=13,m=0)
}

#Baseline Correction
LNCaP.IF <- baseline(filtLNCaP.IF, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
LNCaP.IF <- LNCaP.IF@corrected

ID=as.character(c(vec1,vec2))

#NORMALIZE spectra (vector normalization)
LNCaP.0Gy=sweep(LNCaP.0Gy, MARGIN=1, apply(LNCaP.0Gy,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))
LNCaP.IF=sweep(LNCaP.IF, MARGIN=1, apply(LNCaP.IF,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))

dataall=rbind(LNCaP.0Gy,LNCaP.IF)

#PCA
library("factoextra")
data.pca <- prcomp(dataall ,center = TRUE,scale = FALSE)

fviz_pca_ind(data.pca,label="none",title="",habillage = ID, addEllipses = TRUE, ellipse.level = 0.95)+
  theme_minimal()+
  #xlab("PC1 (52.9%)")+
  #ylab("PC2 (19.3%)")+
  scale_color_brewer(palette="Set1")

##########################################################
#PNT1A

###################First two data sets
data=readMat("PNT1A_0Gy.mat") #This loads the matlab file into the workspace as a 'list' object
PNT1A.0Gy=data$data #Extracting the PNT1A PNT1A.0Gytral data matrix

data=readMat("PNT1A_OF.mat") #This loads the matlab file into the workspace as a 'list' object
PNT1A.OF=data$data 

#SG Filtering 
numrow=nrow(PNT1A.0Gy)
numcol=ncol(PNT1A.0Gy)

filtPNT1A.0Gy=as.matrix(PNT1A.0Gy,nrow=numrow,ncol=numcol)
library(signal)

for (i in 1:numrow){
  filtPNT1A.0Gy[i,]=sgolayfilt(PNT1A.0Gy[i,],p=5,n=13,m=0)
}
#Baseline Correction
PNT1A.0Gy <- baseline(filtPNT1A.0Gy, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
PNT1A.0Gy <- PNT1A.0Gy@corrected

#SG Filtering 
numrow1=nrow(PNT1A.OF)
numcol1=ncol(PNT1A.OF)

filtPNT1A.OF=as.matrix(PNT1A.OF,nrow=numrow1,ncol=numcol1)

for (i in 1:numrow1){
  filtPNT1A.OF[i,]=sgolayfilt(PNT1A.OF[i,],p=5,n=13,m=0)
}
#Baseline Correction
PNT1A.OF <- baseline(filtPNT1A.OF, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
PNT1A.OF <- PNT1A.OF@corrected

#NORMALIZE spectra (vector normalization)
PNT1A.0Gy=sweep(PNT1A.0Gy, MARGIN=1, apply(PNT1A.0Gy,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))
PNT1A.OF=sweep(PNT1A.OF, MARGIN=1, apply(PNT1A.OF,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))

#######################Second two data sets
data=readMat("PNT1A_0Gy_ICCM.mat") #This loads the matlab file into the workspace as a 'list' object
PNT1A.0Gy.ICCM=data$data #Extracting the spectral data matrix

data=readMat("PNT1A_OF_ICCM.mat") #This loads the matlab file into the workspace as a 'list' object
PNT1A.OF.ICCM=data$data 

#SG Filtering
numrow=nrow(PNT1A.0Gy.ICCM)
numcol=ncol(PNT1A.0Gy.ICCM)

filtPNT1A.0Gy.ICCM=as.matrix(PNT1A.0Gy.ICCM,nrow=numrow,ncol=numcol)
library(signal)

for (i in 1:numrow){
  filtPNT1A.0Gy.ICCM[i,]=sgolayfilt(PNT1A.0Gy.ICCM[i,],p=5,n=13,m=0)
}
#Baseline Correction
PNT1A.0Gy.ICCM <- baseline(filtPNT1A.0Gy.ICCM, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
PNT1A.0Gy.ICCM <- PNT1A.0Gy.ICCM@corrected

#SG Filtering
numrow1=nrow(PNT1A.OF.ICCM)
numcol1=ncol(PNT1A.OF.ICCM)

filtPNT1A.OF.ICCM=as.matrix(PNT1A.OF.ICCM,nrow=numrow1,ncol=numcol1)

for (i in 1:numrow1){
  filtPNT1A.OF.ICCM[i,]=sgolayfilt(PNT1A.OF.ICCM[i,],p=5,n=13,m=0)
}
#Baseline Correction
PNT1A.OF.ICCM <- baseline(filtPNT1A.OF.ICCM, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
PNT1A.OF.ICCM <- PNT1A.OF.ICCM@corrected

#NORMALIZE spectra (vector normalization)
PNT1A.0Gy.ICCM=sweep(PNT1A.0Gy.ICCM, MARGIN=1, apply(PNT1A.0Gy.ICCM,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))
PNT1A.OF.ICCM=sweep(PNT1A.OF.ICCM, MARGIN=1, apply(PNT1A.OF.ICCM,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))

##########################################################
#FOUR WAY PCA PNT1A in each of 4 classes 
numrow=nrow(PNT1A.0Gy)
vec1=NULL
vec1[seq(1,numrow)]="PNT1A 0Gy"
numrow2=nrow(PNT1A.0Gy.ICCM)
vec2=NULL
vec2[seq(1,numrow2)]="PNT1A 0Gy + ICCM"
numrow3=nrow(PNT1A.OF)
vec3=NULL
vec3[seq(1,numrow3)]="PNT1A OF"
numrow4=nrow(PNT1A.OF.ICCM)
vec4=NULL
vec4[seq(1,numrow4)]="PNT1A OF + ICCM"

dataall=rbind(PNT1A.0Gy,PNT1A.0Gy.ICCM,PNT1A.OF,PNT1A.OF.ICCM)
ID=as.character(c(vec1,vec2,vec3,vec4))

#PCA 
library("factoextra")
data.pca <- prcomp(dataall ,center = TRUE,scale = FALSE)
variance=data.pca$sdev*data.pca$sdev
pvar=variance/sum(variance)
totvar=sum(pvar[1:5])

fviz_screeplot(data.pca,ncp=10)

fviz_pca_ind(data.pca,label="none",axes=c(1,2),title="",habillage = ID, addEllipses = TRUE, ellipse.level = 0.95)+
  theme_minimal()+
  xlab("PC1 (63.6%)")+
  ylab("PC2 (24.1%)")+
  scale_color_brewer(palette="Set1")

fviz_pca_ind(data.pca,label="none",axes=c(2,3),title="",habillage = ID, addEllipses = TRUE, ellipse.level = 0.95)+
  theme_minimal()+
  #xlab("PC1 (63.6%)")+
  #ylab("PC2 (24.1%)")+
  scale_color_brewer(palette="Set1")

fviz_pca_ind(data.pca,label="none",axes=c(3,4),title="",habillage = ID, addEllipses = TRUE, ellipse.level = 0.95)+
  theme_minimal()+
  #xlab("PC1 (63.6%)")+
  #ylab("PC2 (24.1%)")+
  scale_color_brewer(palette="Set1")


#Plot loadings
loadings = as.data.frame(data.pca$rotation)#Loadings matrix
wavenumber=as.numeric(wavenumber)

#Plot mean spectra
PNT1A_0Gy <- new ("hyperSpec", spc = PNT1A.0Gy, wavelength=wavenumber, data=NULL,label=list("Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))

LDPlot<-ggplot(loadings, aes(x=wavenumber,y=loadings$PC1))
LDPlot+geom_line()+theme_minimal()
LDPlot<-ggplot(loadings, aes(x=wavenumber,y=loadings$PC2))
LDPlot+geom_line()+theme_minimal()
LDPlot<-ggplot(loadings, aes(x=wavenumber,y=loadings$PC3))
LDPlot+geom_line()+theme_minimal()

PNT1A_0Gy.ICCM <- new ("hyperSpec", spc = PNT1A.0Gy.ICCM, wavelength=wavenumber, data=NULL,label=list("Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))
PNT1A_OF <- new ("hyperSpec", spc = PNT1A.OF, wavelength=wavenumber, data=NULL,label=list("Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))
PNT1A_OF.ICCM <- new ("hyperSpec", spc = PNT1A.OF.ICCM, wavelength=wavenumber, data=NULL,label=list("Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))

par(mfrow=c(4,1))
plot(PNT1A_0Gy,"spcmeansd",col=1) #Plotting spectra with mean and sd error bars
plot(PNT1A_0Gy.ICCM,"spcmeansd",col=2) #Plotting spectra with mean and sd error bars
plot(PNT1A_OF,"spcmeansd",col=3) #Plotting spectra with mean and sd error bars
plot(PNT1A_OF.ICCM,"spcmeansd",col=4) #Plotting spectra with mean and sd error bars