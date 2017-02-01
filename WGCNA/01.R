
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\coexpression\\SimulatedData')
library(WGCNA)
options(stringAsFactor=F)

# Here are input parameters of the simulation model
# number of samples or microarrays in the training data
no.obs=50
# now we specify the true measures of eigengene significance
# recall that ESturquoise=cor(y,MEturquoise)
ESturquoise=0; ESbrown= -.6;
ESgreen=.6;ESyellow=0
# Note that we dont specify the eigengene significance of the blue module
# since it is highly correlated with the turquoise module.
ESvector=c(ESturquoise,ESbrown,ESgreen,ESyellow)
# number of genes
nGenes1=3000
# proportion of genes in the turquoise, blue, brown, green, and yellow module #respectively.
simulateProportions1=c(0.2,0.15, 0.08, 0.06, 0.04)
# Note that the proportions dont add up to 1. The remaining genes will be colored grey,
# ie the grey genes are non-module genes.
# set the seed of the random number generator. As a homework exercise change this seed.
set.seed(1)
#Step 1: simulate a module eigengene network.
# Training Data Set I
MEgreen=rnorm(no.obs)
scaledy=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(no.obs)
y=ifelse( scaledy>median(scaledy),2,1)
MEturquoise= ESturquoise*scaledy+sqrt(1-ESturquoise^2)*rnorm(no.obs)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue= .6*MEturquoise+ sqrt(1-.6^2) *rnorm(no.obs)
MEbrown= ESbrown*scaledy+sqrt(1-ESbrown^2)*rnorm(no.obs)
MEyellow= ESyellow*scaledy+sqrt(1-ESyellow^2)*rnorm(no.obs)
ModuleEigengeneNetwork1=data.frame(y,MEturquoise,MEblue,MEbrown,MEgreen, MEyellow)

signif(cor(ModuleEigengeneNetwork1, use="p"),2)

dat1=simulateDatExpr5Modules(MEturquoise=ModuleEigengeneNetwork1$MEturquoise,
                             MEblue=ModuleEigengeneNetwork1$MEblue,
                             MEbrown=ModuleEigengeneNetwork1$MEbrown,
                             MEyellow=ModuleEigengeneNetwork1$MEyellow,
                             MEgreen=ModuleEigengeneNetwork1$MEgreen,
                             nGenes=nGenes1,
                             simulateProportions=simulateProportions1)

datExpr = dat1$datExpr;
truemodule = dat1$truemodule;
datME = dat1$datME;
attach(ModuleEigengeneNetwork1)

table(truemodule)
dim(datExpr)

datExpr=data.frame(datExpr)
ArrayName=paste("Sample",1:dim(datExpr)[[1]], sep="" )
# The following code is useful for outputting the simulated data
GeneName=paste("Gene",1:dim(datExpr)[[2]], sep="" )
dimnames(datExpr)[[1]]=ArrayName
dimnames(datExpr)[[2]]=GeneName

rm(dat1); collectGarbage();
# The following command will save all variables
save.image("Simulated-dataSimulation.RData");

