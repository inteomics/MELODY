########################
##R functions to perform MELODY:
##MEdiation analysis framework in LOgistic regression for high-Dimensional mediators and a binarY outcome;
##based on: Sunyi Chi, Peng Wei and Xuelin Huang (2024) MELODY: Mediation Analysis in Logistic Regression for High-Dimensional Mediators and a Binary Outcome. Manuscript. 
##Created by Sunyi Chi, PhD
##Maintained by Peng Wei, PhD (pwei2@mdanderson.org)

library('SIS')
library('devtools')
library(MASS)
require(msm)
require(survival)
library(glmnet)
library(ncvreg)
library(openxlsx) 
library(DescTools)

#############
##dat is a 5204x17882 dataframe including the bianry outcome variable "chd",
##17873 genes' expression profiles, exposure variable "SEX", 
##and covariates smoking and age: "smok", "AGE1"

sum(dat[,"chd"]) #181
sample_size=dim(dat)[1] #5204
med=dat[,2:17874]

#sex as the primary exposure of interest;
exp=as.matrix(dat[,'SEX'])
colnames(exp)="gender"
Z=scale(cbind(dat[,'smok'],dat[,'AGE1']))
colnames(Z)=c("smoke","age")

Y=dat[,"chd"]
X=exp
M=med
trainrate=0.75
set.seed(111)
train=1:(sample_size*trainrate)
test=(sample_size*trainrate+1):sample_size
######## mediator selection #########
#step 1&2
yt=Y[train]
xt=scale(cbind(X,Z,M)[train,])

model=SIS(xt, yt, family='binomial',penalty = 'MCP',tune='cv',iter=FALSE)
#selected by step 1
#sis.ix0 is index, exposure is 1 cov is 2 and 3 in sis.ix0 so we minus 3 for index to get the mediators' index
SIS_M=model$sis.ix0-3
SIS_M=SIS_M[SIS_M>0]
SIS_M=colnames(M)[SIS_M]
print(SIS_M)

#selected by step 2
xtt=scale(cbind(X[train],Z[train,],M[train, SIS_M]))
fit2=cv.ncvreg(xtt,yt,penalty="MCP", family='binomial')
#The coefficients corresponding to the value of λ that minimizes the cross-validation error can be obtained via coef
ID_p_non <- which(coef(fit2) != 0)
names(ID_p_non)
MCP_M<-names(ID_p_non)[-1]
M_MCP <- M[train, MCP_M]
MCP_M


#step 3
cal_alpha_simple<-function(x){
  data2<-data.frame(Med=x,envir=X[train])
  l<-summary(stats::lm('Med ~.', data=data2))
  invalphap<-(l$coefficients)['envir','Pr(>|t|)']
  return(invalphap)
}

if(length(MCP_M)==1) {inv.alpha.p<-cal_alpha_simple(M_MCP)}else{
  inv.alpha.p<-apply(M_MCP,2, cal_alpha_simple)}
#order <- order(stats::p.adjust(inv.alpha.p, method='fdr'))[1:round(init.cutoff *ncol(Med)) ]
FDR_M=names(which(stats::p.adjust(inv.alpha.p, method='fdr')<0.2))
print(FDR_M)

########## fit model ########
M=scale(med[-train,FDR_M])
Y=Y[-train]
X=scale(X[-train])
cov=scale(Z[-train,])

glm_XM <- glm(Y ~ X+cov+M, family=binomial)


R2_XM=PseudoR2(glm_XM,"McFadden")
coef_b = glm_XM$coefficients[-c(1:4)]
coef_r = glm_XM$coefficients[2]

glm_M <- glm(Y ~ M+cov, family=binomial)
R2_M=PseudoR2(glm_M,"McFadden")
med1 <- lm(M ~ X)
coef_a = med1$coefficients[2,]

glm_X <- glm(Y ~ X+cov, family=binomial)
coef_c = glm_X$coefficients[2]
R2_X=PseudoR2(glm_X,"McFadden")
glm_Z <- glm(Y ~ cov, family=binomial)
R2_Z=PseudoR2(glm_Z,"McFadden")

R2_med  <- (R2_X  + R2_M  - R2_XM- R2_Z) /(1- R2_Z)
SOS=R2_med/(R2_X- R2_Z) /(1- R2_Z)
product=c(coef_a%*%coef_b/coef_c,coef_a%*%coef_b,coef_c-coef_r,(coef_c-coef_r)/coef_c)
names(product)=c("ab/c","ab","c-r","(c-r)/c")
options(digits=3)
(R2_X- R2_Z) /(1- R2_Z)  
(R2_M- R2_Z) /(1- R2_Z)  
(R2_XM- R2_Z) /(1- R2_Z)
R2_Z
R2_med
SOS
product


