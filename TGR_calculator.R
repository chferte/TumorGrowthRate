## Charles Ferté
## Institut Gustave Roussy, Villejuif, France
## Asessment of Tumor Growth Rates TGR) in the clinical setting
## Performed under the supervision of Serge Koscielny, Bernard Escudier and Jean-Charles Soria
## With the great help of Antoine Hollebecque (Gustave Roussy), Mariana Fernandez (Gustave Roussy), Christophe Massard (Gustave Roussy) 
## and Brian Bot (Sage Bionetworks, Seattle)

########################################################################
# definition of TGR
########################################################################
# Tumor size (D) was defined as the sum of the longest diameters of the target lesions 
# as per the Response Evaluation Criteria in Solid Tumors (RECIST) criteria [1]. 
# Let t be the time expressed in months at the tumor evaluation. 
# Assuming the tumor growth follows an exponential law, 
# Vt the tumor volume at time t is equal to Vt=V0 exp(TG.t), 
# where V0 is volume at baseline, and TG is the growth rate. 
# We approximated the tumor volume (V) by V = 4 π R3 / 3, where R, the radius of the sphere is equal to D/2. 
# Consecutively, TG is equal to TG=3 Log(Dt/D0)/t, where D0 and Dt are the diameters (RECIST sums) at baseline and t, respectively.

# To report the tumor growth rate (TGR) results in a clinically meaningful way, 
# we expressed TGR as a percent increase in tumor volume during one month using the following transformation: 
# TGR = 100 (exp(TG) -1), where exp(TG) represents the exponential of TG.

# We calculated the TGR across clinically relevant treatment periods (Figure 1): 
# (i) TGR REFERENCE assessed during the wash-out period (off-therapy) before the introduction of the experimental drug, 
#(ii) TGR EXPERIMENTAL assessed during the first cycle of treatment (i.e.: between the drug introduction and the first evaluation, on-therapy),

########################################################################
#  Notable references about TGR:
########################################################################
# 1: A simulation model of the natural history of human breast cancer. 
# Koscielny S, Tubiana M, Valleron AJ. 
# Br J Cancer. 1985 Oct;52(4):515-24. PubMed PMID: 4063132; PubMed Central PMCID: PMC1977243.

# 2: Tumour growth rates and RECIST criteria in early drug development. Gomez-Roca C, Koscielny S, Ribrag V, Dromain C, Marzouk I, Bidault F, Bahleda R, Ferté C, Massard C, Soria JC. 
# Eur J Cancer. 2011 Nov;47(17):2512-6. doi:10.1016/j.ejca.2011.06.012. Epub 2011 Jul 15. PubMed PMID: 21763126.

# 3: Tumor Growth Rate (TGR) provides useful information to evaluate Sorafenib and Everolimus treatment in metastatic renal cell carcinoma (mRCC) patients. An integrated analysis of the TARGET and RECORD phase III trials data.
# Ferté C, Koscielny S, Albiges L, Rocher L, Soria JC, Iacovelli R, Loriot Y, Fizazi K, Escudier B.
# presented as Posted Discussion, GU Session, ASCO Annual meeting 2012 
# submitted for publication

# 4: Tumor Growth Rate (TGR) provides useful information for patients enrolled in phase I trials 
# and yields clear specific drug profiles.
# Ferté C, et al (manuscript in preparation)
# presented as Posted Discussion, Developmental Therapeutics Session, ASCO Annual meeting 2013 


#########################################################################################################################
# load the data
#########################################################################################################################

# you are invited to input your own file (e.g. MyData.txt)

# requirement is a table excel spreadsheet converted into a tabe delimited .txt file
# must contain the follwoing columns names: "RECIST_BASELINE","RECIST_BEFORE","RECIST_EVAL1","SCANDATE_BASELINE","SCANDATE_BEFORE","SCANDATE_EVAL1"
# with numeric values for the following columns: "RECIST_BASELINE" "RECIST_BEFORE" "RECIST_EVAL1"
# with dates entered as mm/dd/yyyy for the following columns: "SCANDATE_BASELINE" "SCANDATE_BEFORE" "SCANDATE_EVAL1"  


# # for clarity purposes, as example, we point to an example of such .txt file 
# # this file is available through Synapse. 
# # You'll love it: it's super useful and free ! here is some general information about it: http://sagebase.org/synapse-overview/
# # just create a synapse account online @ www.synapse.org 
# # then, install the R client package with the two command lines:
# source("http://depot.sagebase.org/CRAN.R")
# pkgInstall("synapseClient")
# # then run the following lines to get the example.txt file
# library(synapseClient)
# synapseLogin(username="insert here your login",password="insert here your password")
# myFile <- synGet("syn1936427")
# myFile <- myFile@filePath


myData <- read.table(file=myFile)

##########################################################################################################
#  Compute the tumor Growth Rates TGR.ref and TGR.exp (TGR ref and TGR exp the first cycle, respectively)
##########################################################################################################

# first define the reference period and the experimental period (in months)
myData$ref.period <- as.numeric(difftime(myData$SCANDATE_BASELINE,myData$SCANDATE_BEFORE))*12/365.25
myData$exp.period <- as.numeric(difftime(myData$SCANDATE_EVAL1,myData$SCANDATE_BASELINE))*12/365.25

# compute the TGR 
myData$TGR.ref <- 100*(exp(3*log(myData$RECIST_BASELINE/myData$RECIST_BEFORE)/(myData$ref.period))-1)
myData$TGR.exp <- 100*(exp(3*log(myData$RECIST_EVAL1/myData$RECIST_BASELINE)/(myData$exp.period))-1)

############################################################################################################
# compare the TGR ref and the TGR exp (Pairwise comparison by wilcoxon ranked signed test)
# (note that this is a pairwise comparison since each patient is used as his/her own control)
############################################################################################################

# basic descriptive statistics
summary(myData$TGR.ref)
summary(myData$TGR.exp)

# comparison of the two periods
wilcox.test(myData$TGR.ref,myData$TGR.exp,paired=TRUE)

##########################################################################
# plot the TGR ref and the TGR exp 
##########################################################################
par(mfrow=c(1,1))
plot(myData$TGR.ref,myData$TGR.exp,xlim=c(-200,200),ylim=c(-200,200),pch=20,cex=.9,col="gray60",
     axes=FALSE, xlab= "TGR Reference", ylab="TGR Experimental",cex.lab=.8)
axis(1,las=1,at=c(-200,-100,0,100,200),lab=c("-200 %","-100 %","0 %","100 %","200 %"),cex.axis=.7, font=2)
axis(2,las=1,at=c(-200,-100,0,100,200),lab=c("-200 %","-100 %","0 %","100 %","200 %"),cex.axis=.7, font=2)
abline(h=0,v=0)
abline(coef=as.vector(c(0,1)), col="orange", lty=2,lwd=2.5)
text(x=-130, y=-180,labels=paste("Pairwise comparison:\np value =",
    format(wilcox.test(myData$TGR.ref,myData$TGR.exp,paired=TRUE)$p.value,digits=3)),adj=0,cex=.8)
text(-200,-100,"orange line set for: \nTGR ref = TGR exp",col="orange", cex=.55, font=4,adj=0)
text(x=130,y=-105,'DECREASE in TGR\n "Antitumor activity"',cex=.9,font=4,col="darkgreen")
text(x=-90,y=100,'INCREASE in TGR\n "No antitumor activity"',cex=.9,font=4,col="red")
title("Variation of Tumor Growth Rate (TGR)\nacross the Reference and Experimental periods", font=2)

############################################################################################################
# see Pubmed for recent publications using TGR in oncology:
############################################################################################################

# 1: Ferte C, Fernandez M, Hollebecque A, Koscielny S, Levy A, Massard C, Balheda
# R, Bot BM, Gomez Roca C, Dromain C, Ammari S, Soria JC. Tumor Growth Rate (TGR)
# is an early indicator of anti-tumor drug activity in phase I clinical trials.
# Clin Cancer Res. 2013 Nov 22. [Epub ahead of print] PubMed PMID: 24240109.


# 2: Ferté C, Koscielny S, Albiges L, Rocher L, Soria JC, Iacovelli R, Loriot Y,
# Fizazi K, Escudier B. Tumor Growth Rate Provides Useful Information to Evaluate
# Sorafenib and Everolimus Treatment in Metastatic Renal Cell Carcinoma Patients:
# An Integrated Analysis of the TARGET and RECORD Phase 3 Trial Data. Eur Urol.
# 2013 Aug 15. doi:pii: S0302-2838(13)00831-2. 10.1016/j.eururo.2013.08.010. [Epub 
# ahead of print] PubMed PMID: 23993162.

