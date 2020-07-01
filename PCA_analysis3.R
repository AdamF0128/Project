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
##########################################################
#LNCaP
data=readMat("LNCaP_0Gy.mat") #This loads the matlab file into the workspace as a 'list' object
LNCaP.0Gy=data$data #Extracting the LNCaP spectral data matrix

data=readMat("LNCaP_IF.mat") #This loads the matlab file into the workspace as a 'list' object
LNCaP.IF=data$data 

