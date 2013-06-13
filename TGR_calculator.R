## Charles Ferté, MD PhD candidate
## Institut Gustave Roussy, Villejuif, France
## Asessment of Tumor Growth Rates TGR) in the clinical setting
## Performed under the supervision of Serge Koscielny, PhD & Jean-Charles Soria, MD PhD


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
# presented as Posted Discussion, GU Session, ASCO annual meeting 2012 
# submitted for publication

# 4: Tumor Growth Rate (TGR) provides useful information for patients enrolled in phase I trials 
# and yields clear specific drug profiles.
# Ferté C, et al (manuscript in preparation)
# presented as Posted Discussion, Developmental Therapeutics Session, ASCO annual meeting 2012 


########################################################################
# load the data
########################################################################

myFile <- "/volumes/Macintosh HD/Users/chferte/Documents/TRAVAUX/TGR_SITEP/working_database/curated_dataset_TGR_SITEP.txt"
S1 <- read.delim2(file=myFile,na.strings=c("NA", "", " "), as.is=T,header=T)
rownames(S1) <- S1$Patient_Name
rm(myFile)
S1$ti <- as.numeric(S1$ti)
S1$tf <- as.numeric(S1$tf)

# remove the patients that did not receive the treatment (screen failures) (n=76)
table(is.na(S1$Date_Of_C1D1))
S1 <- S1[!is.na(S1$Date_Of_C1D1),]

#remove the patients with non measurable disease at baseline (n=48)
table(is.na(S1$RECIST_BASELINE))
S1 <- S1[!is.na(S1$RECIST_BASELINE),]

S1 <- 
########################################################################
#  Compute the tumor Growth Rates TGR.ref and TGR.exp (TGR ref and TGR exp the first cycle, respectively)
########################################################################

S1$TGR.ref <- 100*(exp(3*log(S1$RECIST_BASELINE/S1$RECIST_BEFORE)/(S1$ti))-1)
S1$TGR.exp <- 100*(exp(3*log(S1$RECIST_EVAL1/S1$RECIST_BASELINE)/(S1$tf))-1)

##########################################################################
# plot the TGR ref and the TGR exp 
##########################################################################
par(mfrow=c(1,1))
plot(TGR.ref,TGR.exp,xlim=c(-200,200),ylim=c(-200,200),pch=20,cex=.9,col="gray60",
     axes=FALSE, xlab= "TGR reference period", ylab="TGR experimental period",cex.lab=.7)
axis(1,las=1,at=c(-200,-100,0,100,200),lab=c("-200 %","-100 %","0 %","100 %","200 %"),cex.axis=.6)
axis(2,las=1,at=c(-200,-100,0,100,200),lab=c("-200 %","-100 %","0 %","100 %","200 %"),cex.axis=.6)
abline(h=0,v=0)
abline(coef=as.vector(c(0,1)), col="orange", lty=2,lwd=2.5)
#text(x=-130, y=-180,labels=paste("Pairwise comparison:\np value =",
#    format(wilcox.test(TGR.ref,TGR.exp,paired=TRUE)$p.value,digits=3)),adj=0,cex=.8)
text(-180,-195,"TGR ref = TGR exp", srt=44,col="orange", cex=.55, font=4,adj=0)
text(x=130,y=-105,'DECREASE in TGR\n "Antitumor activity"',cex=.9,font=4,col="darkgreen")
text(x=-90,y=100,'INCREASE in TGR\n "No antitumor activity"',cex=.9,font=4,col="red")
