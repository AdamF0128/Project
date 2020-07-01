#LOADING THE DATA
rm(list=ls())
dev.off(dev.list()["RStudioGD"]) #Clear all plots

# install.packages("R.matlab")
# install.packages("hyperLNCaP_0Gy")
# install.packages("ChemometricsWithR")
# install.packages("baseline")
# install.packages("gridExtra")
library(R.matlab)   #These lines load the library files into the workspace
library(hyperSpec)
library(ChemometricsWithR)
library(baseline)
library(ggplot2)
library(gridExtra)
library(signal)

data=readMat("LNCaP_0Gy.mat") #This loads the matlab file into the workspace as a 'list' object
LNCaP_0Gy=data$data #Extracting the spectral data matrix

#Loading wavenumber data
wavenumber=readMat("wavenumber.mat")
wavenumber=wavenumber$wavenumber

#SG Filtering
numrow=nrow(LNCaP_0Gy)
numcol=ncol(LNCaP_0Gy)
filtLNCaP_0Gy=as.matrix(LNCaP_0Gy,nrow=numrow,ncol=numcol)

for (i in 1:numrow){
  filtLNCaP_0Gy[i,]=sgolayfilt(LNCaP_0Gy[i,],p=5,n=13,m=0)
}

LNCaP.0Gy=sweep(filtLNCaP_0Gy, MARGIN=1, apply(filtLNCaP_0Gy,1,function(x) sqrt(sum(x^2))) , FUN="/") #Vector normalisation function(x)sqrt(sum(x^2))

LNCaP_0Gy.als <- baseline(LNCaP.0Gy, lambda = 9.5, p = 0.2, maxit = 20, method='als')
LNCaP_0Gy.als <- LNCaP_0Gy.als@corrected

LNCaP_0Gy.LP <- baseline(LNCaP.0Gy, lambda=2, hwi=210,  it=5,  int=1050,method='fillPeaks')
LNCaP_0Gy.LP <- LNCaP_0Gy.LP@corrected

LNCaP_0Gy.MP <- baseline(LNCaP.0Gy, degree=4, method='modpolyfit')
LNCaP_0Gy.MP <- LNCaP_0Gy.MP@corrected

wavenumber=as.numeric(wavenumber)

par(mfrow=c(3,1))

LNCaP_0Gy_als <- new ("hyperSpec", spc = LNCaP_0Gy.als, wavelength=wavenumber, data=NULL,label=list(spc="Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))
plot(LNCaP_0Gy_als,"spcmeansd",col=3) #Plotting spectra with mean and sd error bars

LNCaP_0Gy_LP <- new ("hyperSpec", spc = LNCaP_0Gy.LP, wavelength=wavenumber, data=NULL,label=list(spc="Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))
plot(LNCaP_0Gy_LP,"spcmeansd",col=4) #Plotting spectra with mean and sd error bars

LNCaP_0Gy_MP <- new ("hyperSpec", spc = LNCaP_0Gy.MP, wavelength=wavenumber, data=NULL,label=list(spc="Intensity (a.u.)",.wavelength=expression(Delta * tilde(nu)/cm^-1)))
plot(LNCaP_0Gy_MP,"spcmeansd",col=5) #Plotting spectra with mean and sd error bars
