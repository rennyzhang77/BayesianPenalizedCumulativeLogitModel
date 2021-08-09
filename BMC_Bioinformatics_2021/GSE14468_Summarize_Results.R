rm(list=ls())
gc()
options(expressions=500000)
setwd("your-directory-that-stores-GSE14468-model-data_objects.RData")
library(Biobase)

####### MODEL I ##################
load("Model1BMC.RData")
end.time.1-start.time.1
### select variable using bayes factor based on gammaXbeta
summary(BetgmaFactor)
### select variable using 95% equal-tailed credible interval of gammaXbeta
sum(e.hd.1)
### select variable using 95% HPD credible interval of gammaXbeta
sum(e.ci.1)
rm(list=ls())

####### MODEL II ##################
## Summarize Model 2
## MODEL 2 Prior fixed at pi=0.50
load("Model2FixedPi50BMC.RData")
end.time.2.pi.50-start.time.2.pi.50
length(gamma.index.2)
length(BFbeta.2)
model.2.gamma.pi.beta.5<-featureNames(expression)[gamma.index.2]
model.2.gamma.pi.beta.5
model.2.beta.pi.beta.5<-featureNames(expression)[BFbeta.2]
sum(gamma.index.2%in%BFbeta.2)
sum(e.ci.2)
sum(e.hd.2) 
rm(list=ls()[-grep("model.",ls())])

## MODEL 2 Prior fixed at pi=0.05
load("Model2FixedPi05BMC.RData")
end.time.2.pi.05-start.time.2.pi.05
length(gamma.index.2)
length(BFbeta.2)
model.2.gamma.pi.beta.05<-featureNames(expression)[gamma.index.2]
model.2.gamma.pi.beta.05
model.2.beta.pi.beta.05<-featureNames(expression)[BFbeta.2]
sum(gamma.index.2%in%BFbeta.2)
sum(e.ci.2)
sum(e.hd.2) 
rm(list=ls()[-grep("model.",ls())])

## Model 2 Pi ~ Beta(1,19)
load("Model2Beta1_19BMC.RData")
end.time.2.beta.1.19-start.time.2.beta.1.19
length(gamma.index.2)
length(BFbeta.2)
model.2.gamma.pi.beta.1.19<-featureNames(expression)[gamma.index.2]
model.2.gamma.pi.beta.1.19
model.2.beta.pi.beta.1.19<-featureNames(expression)[BFbeta.2]
sum(gamma.index.2%in%BFbeta.2)
sum(e.ci.2)
sum(e.hd.2) 
rm(list=ls()[-grep("model.",ls())])

## Extra Model II Prior is Beta(0.01, 0.19)
load("Model2Beta001_019BMC.RData")
end.time.2.beta.001.019-start.time.2.beta.001.019
length(gamma.index.2)
length(BFbeta.2)
model.2.gamma.pi.beta.001.019<-featureNames(expression)[gamma.index.2]
model.2.gamma.pi.beta.001.019
model.2.beta.pi.beta.001.019<-featureNames(expression)[BFbeta.2]
sum(gamma.index.2%in%BFbeta.2)
sum(e.ci.2)
sum(e.hd.2) 
rm(list=ls()[-grep("model.",ls())])

####### MODEL III ##################
### Summarize MODEL 3
## Model 3 Pi fixed at 0.50
load("Model3FixedPi50BMC.RData")
end.time.3.pi.50-start.time.3.pi.50
length(gamma.index.3)
length(BFbeta.3)
model.3.gamma.pi.5<-featureNames(expression)[gamma.index.3]
model.3.gamma.pi.5
model.3.beta.pi.5<-featureNames(expression)[BFbeta.3]
sum(gamma.index.3%in%BFbeta.3)
sum(e.ci.3)
sum(e.hd.3) 
rm(list=ls()[-grep("model.",ls())])

## Model 3 Pi fixed at 0.05
load("Model3FixedPi05BMC.RData")
end.time.3.pi.05-start.time.3.pi.05
length(gamma.index.3)
length(BFbeta.3)
model.3.gamma.pi.05<-featureNames(expression)[gamma.index.3]
model.3.gamma.pi.05
model.3.beta.pi.05<-featureNames(expression)[BFbeta.3]
sum(gamma.index.4%in%BFbeta.3)
sum(e.ci.3)
sum(e.hd.3) 
rm(list=ls()[-grep("model.",ls())])

## Model 3 Pi ~ Beta(1,19)
load("Model3Beta1_19BMC.RData")
end.time.3.beta.1.19 - start.time.3.beta.1.19
length(gamma.index.3)
length(BFbeta.3)
model.3.gamma.pi.beta.1.19<-featureNames(expression)[gamma.index.3]
model.3.gamma.pi.beta.1.19
model.3.beta.pi.beta.1.19<-featureNames(expression)[BFbeta.3]
sum(gamma.index.3%in%BFbeta.3)
sum(e.ci.3)
sum(e.hd.3)
rm(list=ls()[-grep("model.",ls())])

## Extra Model 3 Prior is Beta(0.01, 0.19)
load("Model3Beta001_019BMC.RData")
end.time.3.beta.001.019 - start.time.3.beta.001.019
length(gamma.index.3)
length(BFbeta.3)
model.3.gamma.pi.beta.001.019<-featureNames(expression)[gamma.index.3]
model.3.gamma.pi.beta.001.019
model.3.beta.pi.beta.001.019<-featureNames(expression)[BFbeta.3]
sum(gamma.index.3%in%BFbeta.3)
sum(e.ci.3)
sum(e.hd.3)
rm(list=ls()[-grep("model.",ls())])

####### MODEL IV ##################
####### Summarize MODEL IV *
## Pi fixed at 0.50
load("Model4FixedPi50BMC.RData")
end.time.4.pi.50-start.time.4.pi.50
length(gamma.index.4)
length(BFbeta.4)
model.4.gamma.pi.5<-featureNames(expression)[gamma.index.4]
model.4.gamma.pi.5
model.4.beta.pi.5<-featureNames(expression)[BFbeta.4]
sum(gamma.index.4%in%BFbeta.4)
sum(e.ci.4)
sum(e.hd.4) 
rm(list=ls()[-grep("model.",ls())])

load("Model4FixedPi05BMC.RData")
end.time.4.pi.05-start.time.4.pi.05
length(gamma.index.4)
length(BFbeta.4)
model.4.gamma.pi.05<-featureNames(expression)[gamma.index.4]
model.4.gamma.pi.05
model.4.beta.pi.05<-featureNames(expression)[BFbeta.4]
sum(gamma.index.4%in%BFbeta.4)
sum(e.ci.4)
sum(e.hd.4) 
rm(list=ls()[-grep("model.",ls())])

load("Model4Beta1_19BMC.RData")
end.time.4.beta.1.19-start.time.4.beta.1.19
length(gamma.index.4)
length(BFbeta.4)
model.4.gamma.pi.beta.1.19<-featureNames(expression)[gamma.index.4]
model.4.gamma.pi.beta.1.19
model.4.beta.pi.beta.1.19<-featureNames(expression)[BFbeta.4]
sum(gamma.index.4%in%BFbeta.4)
sum(e.ci.4)
sum(e.hd.4)
rm(list=ls()[-grep("model.",ls())])

## Extra Model IV Prior is Beta(0.01, 0.19)
load("Model4Beta001_019BMC.RData")
end.time.4.beta.001.019-start.time.4.beta.001.019
length(gamma.index.4)
length(BFbeta.4)
model.4.gamma.pi.beta.001.019<-featureNames(expression)[gamma.index.4]
model.4.gamma.pi.beta.001.019
model.4.beta.pi.beta.001.019<-featureNames(expression)[BFbeta.4]
sum(gamma.index.4%in%BFbeta.4)
sum(e.ci.4)
sum(e.hd.4)
rm(list=ls()[-grep("model.",ls())])


##### Boxplots in Supplemental Figures #########################################################
plot.frame<-data.frame(t(exprs(expression)))
cyto.risk<-gsub("cytogenetic ","",pData(rma.cyto.risk)$"risk:ch1")
boxplot(X242520_s_at~cyto.risk, data=plot.frame, xlab="", ylab="ARMH1")
dev2bitmap("ARMH1.pdf","pdfwrite",res=300)
boxplot(X206135_at~cyto.risk, data=plot.frame, xlab="", ylab="ST18")
dev2bitmap("ST18.pdf","pdfwrite",res=300)
boxplot(X208886_at~cyto.risk, data=plot.frame, xlab="", ylab="H1-0")
dev2bitmap("H1-0.pdf","pdfwrite",res=300)
boxplot(X205382_s_at~cyto.risk, data=plot.frame, xlab="", ylab="CFD")
dev2bitmap("CFD.pdf","pdfwrite",res=300)
boxplot(X205844_at~cyto.risk, data=plot.frame, xlab="", ylab="VNN1")
dev2bitmap("VNN1.pdf","pdfwrite",res=300)
boxplot(X1553808_a_at~cyto.risk, data=plot.frame, xlab="", ylab="NKX2-3")
dev2bitmap("NKX2-3.pdf","pdfwrite",res=300)
boxplot(X224596_at~cyto.risk, data=plot.frame, xlab="", ylab="SLC44A1")
dev2bitmap("SLC44A1.pdf","pdfwrite",res=300)
boxplot(X225782_at~cyto.risk, data=plot.frame, xlab="", ylab="MSRB3")
dev2bitmap("MSRB3.pdf","pdfwrite",res=300)
boxplot(X228170_at~cyto.risk, data=plot.frame, xlab="", ylab="OLIG1")
dev2bitmap("OLIG1.pdf","pdfwrite",res=300)
boxplot(X210146_x_at~cyto.risk, data=plot.frame, xlab="", ylab="LILR2")
dev2bitmap("LILR2.pdf","pdfwrite",res=300)
boxplot(X208892_s_at~cyto.risk, data=plot.frame, xlab="", ylab="DUSP6")
dev2bitmap("DUSP6.pdf","pdfwrite",res=300)
boxplot(X201427_s_at~cyto.risk, data=plot.frame, xlab="", ylab="SELENOP")
dev2bitmap("SELENOP.pdf","pdfwrite",res=300)
boxplot(X210997_at~cyto.risk, data=plot.frame, xlab="", ylab="HGF")
dev2bitmap("HGF.pdf","pdfwrite",res=300)
boxplot(X203973_s_at~cyto.risk, data=plot.frame, xlab="", ylab="CEBPD")
dev2bitmap("CEBPD.pdf","pdfwrite",res=300)
boxplot(X204082_at~cyto.risk, data=plot.frame, xlab="", ylab="PBX3")
dev2bitmap("PBX3.pdf","pdfwrite",res=300)
boxplot(X210755_at~cyto.risk, data=plot.frame, xlab="", ylab="HGF")
dev2bitmap("HGF.pdf","pdfwrite",res=300)
boxplot(X210933_s_at~cyto.risk, data=plot.frame, xlab="", ylab="FSCN1")
dev2bitmap("FSCN1.pdf","pdfwrite",res=300)
boxplot(X214329_x_at~cyto.risk, data=plot.frame, xlab="", ylab="TNFSF10")
dev2bitmap("TNFSF10.pdf","pdfwrite",res=300)
boxplot(X228827_at~cyto.risk, data=plot.frame, xlab="", ylab="RUNX1T1")
dev2bitmap("RUNX1T1.pdf","pdfwrite",res=300)
boxplot(X228904_at~cyto.risk, data=plot.frame, xlab="", ylab="HOXB3")
dev2bitmap("HOXB3.pdf","pdfwrite",res=300)
boxplot(X206283_s_at~cyto.risk, data=plot.frame, xlab="", ylab="TAL1")
dev2bitmap("TAL1.pdf","pdfwrite",res=300)
boxplot(X223044_at~cyto.risk, data=plot.frame, xlab="", ylab="SLC40A1")
dev2bitmap("SLC40A1.pdf","pdfwrite",res=300)
boxplot(X210783_x_at~cyto.risk, data=plot.frame, xlab="", ylab="CLEC11A")
dev2bitmap("CLEC11A.pdf","pdfwrite",res=300)
boxplot(X214651_s_at~cyto.risk, data=plot.frame, xlab="", ylab="214651_s_at")
boxplot(X204961_s_at~cyto.risk, data=plot.frame, xlab="", ylab="204961_s_at")


