## Charles Ferté
## Institut Gustave Roussy, Villejuif, France
## Asessment of Tumor Growth Rates TGR) in the clinical setting
## Performed under the supervision of Serge Koscielny, Bernard Escudier and Jean-Charles Soria
## With the help of Antoine Hollebecque, Mariana Fernandez and Brian Bot

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
# Consecutively, TG is equal to TG=3 Log(Dt/D0)/t. 
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
# requirement is a table excel spreadsheet converted into a tabe delimited .txt file
# with numeric values for the following columns: "RECIST_BASELINE" "RECIST_BEFORE" "RECIST_EVAL1"
# and with dayes entered as mm/dd/yyyy for the following columns: "SCANDATE_BASELINE" "SCANDATE_BEFORE" "SCANDATE_EVAL1"  

########################################################################
myFileURL <- "https://raw.github.com/chferte/TumorGrowthRate/master/TGR_example.txt"
S1 <- read.table(file=myFileURL)
S1 <- read.table("/Users/chferte/Documents/TRAVAUX/TGR_SITEP/TumorGrowthRate/TGR_example.txt")

########################################################################
#  Compute the tumor Growth Rates TGR.ref and TGR.exp (TGR ref and TGR exp the first cycle, respectively)
########################################################################

# first define the reference period and the experimental period (in months)
S1$ref.period <- as.numeric(difftime(S1$SCANDATE_BASELINE,S1$SCANDATE_BEFORE))*12/365.25
S1$exp.period <- as.numeric(difftime(S1$SCANDATE_EVAL1,S1$SCANDATE_BASELINE))*12/365.25

# compute the TGR 
S1$TGR.ref <- 100*(exp(3*log(S1$RECIST_BASELINE/S1$RECIST_BEFORE)/(S1$ref.period))-1)
S1$TGR.exp <- 100*(exp(3*log(S1$RECIST_EVAL1/S1$RECIST_BASELINE)/(S1$exp.period))-1)

##########################################################################
# compaire the TGR ref and the TGR exp (Pairwise comparison by wilcoxon ranked signed test)
# (note that this is a pairwise comparison since each patient is used as his/her own control)
##########################################################################
summary(S1$TGR.ref)
summary(S1$TGR.exp)
wilcox.test(S1$TGR.ref,S1$TGR.exp,paired=TRUE)

##########################################################################
# plot the TGR ref and the TGR exp 
##########################################################################
par(mfrow=c(1,1))
plot(S1$TGR.ref,S1$TGR.exp,xlim=c(-200,200),ylim=c(-200,200),pch=20,cex=.9,col="gray60",
     axes=FALSE, xlab= "TGR reference period", ylab="TGR experimental period",cex.lab=.7)
axis(1,las=1,at=c(-200,-100,0,100,200),lab=c("-200 %","-100 %","0 %","100 %","200 %"),cex.axis=.6)
axis(2,las=1,at=c(-200,-100,0,100,200),lab=c("-200 %","-100 %","0 %","100 %","200 %"),cex.axis=.6)
abline(h=0,v=0)
abline(coef=as.vector(c(0,1)), col="orange", lty=2,lwd=2.5)
text(x=-130, y=-180,labels=paste("Pairwise comparison:\np value =",
    format(wilcox.test(S1$TGR.ref,S1$TGR.exp,paired=TRUE)$p.value,digits=3)),adj=0,cex=.8)
text(-200,-100,"orange line set for: \nTGR ref = TGR exp",col="orange", cex=.55, font=4,adj=0)
text(x=130,y=-105,'DECREASE in TGR\n "Antitumor activity"',cex=.9,font=4,col="darkgreen")
text(x=-90,y=100,'INCREASE in TGR\n "No antitumor activity"',cex=.9,font=4,col="red")
title("Variation of TGR across reference and experimental periods", font=2)


