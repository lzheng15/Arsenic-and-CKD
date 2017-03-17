################################################################################

# ARSENIC & Incident CKD BY ARSENIC METABOLITES
# Original code from: Katherine Moon (Arsenic and CVD) at kmoon9@jhu.edu
# Side Note: Kat got this code from Maria Tellez-Plaza (Cd and CVD)
# Analyst: Laura Zheng 

################################################################################
#                         TABLE OF CONTENTS
# ------------------------------------------------------------------------------

#  1. Setup directories, load data
#  2. Create restricted quadratic splines for all arsenic species and metabolites
#  3. Histograms to find breaks and distributions, and related things 
#  4. Run Cox regressions for incident CKD by arsenic metabolites
#  5. Plot for Arsenic metabolites and Incident CKD  							
#
#		5-A. Plot for Sum of inorganic and methylated species
#		5-B. Plot for inorganic arsenic (Arsenate and Arsenite)
#		5-C. Plot for MMA
#  	5-D. Plot for DMA
#
#
################################################################################

#  1. Setup directories, load data

#call up library to read stata data, call up library for survival data

library(foreign)
library(survival)

#set working directory  in PC and read data into R form, then save, then we can directly reload the R data file
setwd("C:/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr")

#cox data, load in from stata
shs3119people<-read.dta("/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr/shs1-3121-people-cox-130207.dta")
names(shs3119people) <- casefold(names(shs3119people))
save(shs3119people,file="shs3119people.rda")

# read in csv file of data for prospective because the other thing wont work
shs3119people <-read.csv("/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr/130506-test.csv")
save(shs3119people,file="shs3119people.rda")

#load incident survival data
load("C:/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr/shs3119people.rda")

#macbook set working directory

setwd("/Users/LauraZheng/Documents/Hopkins Classes Stuff/Strong Heart Study/basic prospective arsenic gfr")
getwd()

#load on mac 
load("/Users/LauraZheng/Documents/Hopkins Classes Stuff/Strong Heart Study/basic prospective arsenic gfr/shs3119people.rda")
#################################################################################################################################

#  2. Create restricted quadratic splines for all arsenic species and metabolites

#get quantiles first
quantile(shs3119people$sum_over_creat, probs = c(10, 50, 90)/100)
#calculate the splines for sum of inorganic and methylated species
shs3119people$sum.cr.rqs2 <- pmax(log(shs3119people$sum_over_creat)-log(3.8600),0)^2 - pmax(log(shs3119people$sum_over_creat)- log(23.8700),0)^2 #10th, 90th percentiles
shs3119people$sum.cr.rqs3 <- pmax(log(shs3119people$sum_over_creat)-log(9.7200),0)^2 - pmax(log(shs3119people$sum_over_creat)- log(23.8700),0)^2 #50th, 90th percentiles

#to get values of splines we use the pmax function for INORGANIC ARSENIC
#get quantiles first
quantile(shs3119people$asvasiii_ox_cr, probs = c(10, 50, 90)/100)
#calculate the splines for inorganic arsenic 
shs3119people$ias.cr.rqs2 <- pmax(log(shs3119people$asvasiii_ox_cr)-log(0.200),0)^2 - pmax(log(shs3119people$asvasiii_ox_cr)- log(2.612),0)^2 #10th, 90th percentiles
shs3119people$ias.cr.rqs3 <- pmax(log(shs3119people$asvasiii_ox_cr)-log(0.780),0)^2 - pmax(log(shs3119people$asvasiii_ox_cr)- log(2.612),0)^2 #50th, 90th percentiles


#to get values of splines we use the pmax function for MMA
#get quantiles first
quantile(shs3119people$ma_ox_cr, probs = c(10, 50, 90)/100)
#calculate the splines for MMA
shs3119people$mma.cr.rqs2 <- pmax(log(shs3119people$ma_ox_cr)-log(0.48),0)^2 - pmax(log(shs3119people$ma_ox_cr)- log(3.52),0)^2 #10th, 90th percentiles
shs3119people$mma.cr.rqs3 <- pmax(log(shs3119people$ma_ox_cr)-log(1.30),0)^2 - pmax(log(shs3119people$ma_ox_cr)- log(3.52),0)^2 #50th, 90th percentiles

#to get values of splines we use the pmax function for DMA
#get quantiles first
quantile(shs3119people$dma_ox_cr, probs = c(10, 50, 90)/100)
#calculate the splines for DMA
shs3119people$dma.cr.rqs2 <- pmax(log(shs3119people$dma_ox_cr)-log(3.01),0)^2 - pmax(log(shs3119people$dma_ox_cr)- log(18.30),0)^2 #10th, 90th percentiles
shs3119people$dma.cr.rqs3 <- pmax(log(shs3119people$dma_ox_cr)-log(7.31),0)^2 - pmax(log(shs3119people$dma_ox_cr)- log(18.30),0)^2 #50th, 90th percentiles


#########################################################################################################################################

#  3. Histograms to find breaks and distributions, and related things 

# histogram of log of sum of inorganic and methylated as
dist <- hist(log(shs3119people$sum_over_creat))
sumHR.breaks <- dist$breaks      			 
sumHR.breaks <- sumHR.breaks[1:8]					
sumHR.freq <- dist$counts/sum(dist$counts)		 
sumHR.freq <- sumHR.freq[1:7]						
exp(min(sumHR.breaks))
exp(max(sumHR.breaks))

# histogram of log of inorganic arsenic AsV and AsIII
dist.ias<- hist(log(shs3119people$asvasiii_ox_cr))
ias.breaks <- dist.ias$breaks
ias.breaks <- ias.breaks[2:12]
ias.freq <- dist.ias$counts/sum(dist.ias$counts)  
ias.freq <- ias.freq[2:11]
exp(min(ias.breaks))
exp(max(ias.breaks))


# histogram of log of MMA
dist.mma<- hist(log(shs3119people$ma_ox_cr))
mma.breaks <- dist.mma$breaks
mma.breaks <- mma.breaks[3:11]
mma.freq <- dist.mma$counts/sum(dist.mma$counts)  
mma.freq <- mma.freq[3:10]
exp(min(mma.breaks))
exp(max(mma.breaks))

# histogram of log of DMA
dist.dma<- hist(log(shs3119people$dma_ox_cr))
dma.breaks <- dist.dma$breaks
dma.breaks <- dma.breaks[2:9]
dma.freq <- dist.dma$counts/sum(dist.dma$counts)  
dma.freq <- dma.freq[2:8]
exp(min(dma.breaks))
exp(max(dma.breaks))



uas <- seq(log(1.648721),log(54.59815),length=100)  		#look at histogram of log(As), and see where majority of points lie (use this min and max of breaks here) , middle of distribution

ias <- seq(log(0.04978707),log(7.389056),length=100)    	#look at histogram of log(iAs), and see where majority of points lie (use this min and max of breaks here) , middle of distribution

mma <- seq(log(0.1353353),log(7.389056),length=100)      #look at histogram of log(MMA), and see where majority of points lie (use this min and max of breaks here) , middle of distribution

dma <- seq(log(1.648721),log(54.59815),length=100)      #look at histogram of log(DMA), and see where majority of points lie (use this min and max of breaks here) , middle of distribution



######################################################################################

#  4. Run Cox regressions for incident CKD by arsenic metabolites

#sum of inorganic and methylated 
mat1 <- cbind(uas,pmax(uas-(log(3.86)),0)^2-pmax(uas-log(23.87),0)^2,pmax(uas-(log(9.72)),0)^2-pmax(uas-log(23.87),0)^2)  #this is a concatenation thingie
mat1 <- sweep(mat1,2,mat1[25,])     				# "Sweep" out (remove, set equal to zero) predictor at the reference value (10th percentile, log(3.860 ug/g), which is the 25th element in UCD

#IAs
mat2 <- cbind(ias,pmax(ias-(log(0.200)),0)^2-pmax(ias-log(2.612),0)^2,pmax(ias-(log(0.780)),0)^2-pmax(ias-log(2.612),0)^2)  #this is a concatenation thingie
mat2 <- sweep(mat2,2,mat2[28,])         		# "Sweep" out (remove, set equal to zero) predictor at the reference value (10th percentile, log(0.200 ug/g), which is the 25th element in UCD

#MMA 
mat3 <- cbind(mma,pmax(mma-(log(0.48)),0)^2-pmax(mma-log(3.52),0)^2,pmax(mma-(log(1.30)),0)^2-pmax(mma-log(3.52),0)^2)  #this is a concatenation thingie
mat3 <- sweep(mat3,2,mat3[32,])       			# "Sweep" out (remove, set equal to zero) predictor at the reference value (10th percentile, log(0.48 ug/g), which is the 25th element in UCD

#DMA
mat4 <- cbind(dma,pmax(dma-(log(3.01)),0)^2-pmax(dma-log(23.87),0)^2,pmax(dma-(log(7.31)),0)^2-pmax(dma-log(18.30),0)^2)  #this is a concatenation thingie
mat4 <- sweep(mat4,2,mat4[18,])         		# "Sweep" out (remove, set equal to zero) predictor at the reference value (10th percentile, log(3.01 ug/g), which is the 25th element in UCD


#Cox model with splines, for sum of inorganic and methylated
fit1 <- coxph(Surv(entryage,exitage, event=incickd)~sum_ias_cr_ln+sum.cr.rqs2+sum.cr.rqs3+gender+s1bmi+s1edu+ as.factor(s1smoke_n)+s1dmwhb+s1g0+s1htnrx_n+s1sbp+s1gfr+strata(as.factor(location)),data=shs3119people)
summary(fit1)

b3.sum <- summary(fit1)$coeff[1:3,1]   #stores all 3 As coefficients in vector
varb3.sum <- fit1$var[1:3,1:3]  		#stores the variance matrix of the coefﬁcients, for all As coeff.
lnHR.sum <- mat1%*%b3.sum
se.HR.sum <- sqrt(diag(mat1%*%varb3.sum%*%t(mat1))) 
ll.HR.sum <- lnHR.sum-1.96*se.HR.sum
ul.HR.sum <- lnHR.sum+1.96*se.HR.sum


#inorganic arsenic (AsIII and AsV)
fit2 <- coxph(Surv(entryage,exitage, event=incickd)~asvasiii_ox_cr_ln+ias.cr.rqs2+ias.cr.rqs3+gender+s1bmi+s1edu+ as.factor(s1smoke_n)+s1dmwhb+s1g0+s1htnrx_n+s1sbp+s1gfr+strata(as.factor(location)),data=shs3119people)
summary(fit2)

b3.ias <- summary(fit2)$coeff[1:3,1]   #stores all 3 As coefficients in vector
varb3.ias <- fit2$var[1:3,1:3]  		#stores the variance matrix of the coefﬁcients, for all As coeff.
lnHR.ias <- mat2%*%b3.ias
se.HR.ias <- sqrt(diag(mat2%*%varb3.ias%*%t(mat2))) 
ll.HR.ias <- lnHR.ias-1.96*se.HR.ias
ul.HR.ias <- lnHR.ias+1.96*se.HR.ias

#MMA
fit3 <- coxph(Surv(entryage,exitage, event=incickd)~ma_ox_cr_ln+mma.cr.rqs2+mma.cr.rqs3+gender+s1bmi+s1edu+ as.factor(s1smoke_n)+s1dmwhb+s1g0+s1htnrx_n+s1sbp+s1gfr+strata(as.factor(location)),data=shs3119people)
summary(fit3)

b3.mma <- summary(fit3)$coeff[1:3,1]   #stores all 3 As coefficients in vector
varb3.mma <- fit3$var[1:3,1:3]    	#stores the variance matrix of the coefﬁcients, for all As coeff.
lnHR.mma <- mat3%*%b3.mma
se.HR.mma <- sqrt(diag(mat3%*%varb3.mma%*%t(mat3))) 
ll.HR.mma <- lnHR.mma-1.96*se.HR.mma
ul.HR.mma <- lnHR.mma+1.96*se.HR.mma

#MMA + iAs adjustment
fit4 <- coxph(Surv(entryage,exitage, event=incickd)~ma_ox_cr_ln+mma.cr.rqs2+mma.cr.rqs3+gender+s1bmi+s1edu+ as.factor(s1smoke_n)+s1dmwhb+s1g0+s1htnrx_n+s1sbp+s1gfr+asvasiii_ox_cr_ln+strata(as.factor(location)),data=shs3119people)
summary(fit4)


b3.mma2 <- summary(fit4)$coeff[1:3,1]   #stores all 3 As coefficients in vector
varb3.mma2 <- fit4$var[1:3,1:3]      #stores the variance matrix of the coefﬁcients, for all As coeff.
lnHR.mma2 <- mat2%*%b3.mma2
se.HR.mma2 <- sqrt(diag(mat1%*%varb3.mma2%*%t(mat1))) 
ll.HR.mma2 <- lnHR.mma2-1.96*se.HR.mma2
ul.HR.mma2 <- lnHR.mma2+1.96*se.HR.mma2

#DMA
fit5 <- coxph(Surv(entryage,exitage, event=incickd)~dma_ox_cr_ln+dma.cr.rqs2+dma.cr.rqs3+gender+s1bmi+s1edu+ as.factor(s1smoke_n)+s1dmwhb+s1g0+s1htnrx_n+s1sbp+s1gfr+strata(as.factor(location)),data=shs3119people)
summary(fit5)


b3.dma <- summary(fit5)$coeff[1:3,1]   #stores all 3 As coefficients in vector
varb3.dma <- fit5$var[1:3,1:3]      #stores the variance matrix of the coefﬁcients, for all As coeff.
lnHR.dma <- mat2%*%b3.dma
se.HR.dma <- sqrt(diag(mat1%*%varb3.dma%*%t(mat1))) 
ll.HR.dma <- lnHR.dma-1.96*se.HR.dma
ul.HR.dma <- lnHR.dma+1.96*se.HR.dma

#DMA + iAs adjustment
fit6 <- coxph(Surv(entryage,exitage, event=incickd)~dma_ox_cr_ln+dma.cr.rqs2+dma.cr.rqs3+gender+s1bmi+s1edu+ as.factor(s1smoke_n)+s1dmwhb+s1g0+s1htnrx_n+s1sbp+s1gfr+asvasiii_ox_cr_ln+strata(as.factor(location)),data=shs3119people)
summary(fit6)

b3.dma2 <- summary(fit6)$coeff[1:3,1]   #stores all 3 As coefficients in vector
varb3.dma2 <- fit6$var[1:3,1:3]      #stores the variance matrix of the coefﬁcients, for all As coeff.
lnHR.dma2 <- mat2%*%b3.dma2
se.HR.dma2 <- sqrt(diag(mat1%*%varb3.dma2%*%t(mat1))) 
ll.HR.dma2 <- lnHR.dma2-1.96*se.HR.dma2
ul.HR.dma2 <- lnHR.dma2+1.96*se.HR.dma2


####################################################################
##################################################################################
#  5. Plot for Arsenic metabolites and Incident CKD 

# mfrow=c(# rows, # columns) this combines multiple plots into a single graphing window

par(mfrow=c(2,2),mar=c(5,5,5,5))

#(ORDERED FROM LEFT TO RIGHT, TOP TO BOTTOM)  
# 1. First row, first column
# 2. First row, second column
# 3. Second row, first column
# 4. Second row, second column


##############################################
#  	5-A. Plot for Sum of inorganic and methylated species
#par(mar=c(5,5,5,5))
plot(uas,lnHR.sum,type="n",xlim=c(min(sumHR.breaks),max(sumHR.breaks)),ylim=c(log(0.6),log(3.00)),xlab="",ylab="",axes=F)  	#note that type n is blank

rect(sumHR.breaks[1:(length(sumHR.breaks)-1)],rep(log(0.6),length(sumHR.breaks)-1),				# rect(xleft,ybottom,xright,ytop)
     sumHR.breaks[2:length(sumHR.breaks)],log(0.6)+(100/30)*sumHR.freq*(log(3)-log(0.6)),
     density=-1,border="white",col=grey(0.90),lty=1,lwd=1.5)

lines(uas,lnHR.sum, lty=1,lwd=2, col="black")	
lines(uas,ll.HR.sum,lty=2,lwd=2, col="black")	
lines(uas,ul.HR.sum,lty=2,lwd=2, col="black")			

abline(h=0,lty=1,lwd=1) 

axis(1,at=seq(min(sumHR.breaks),max(sumHR.breaks),length=8),labels=round(exp(seq(min(sumHR.breaks),max(sumHR.breaks),length=8)),1),
     cex.axis=1.05,font.axis=1,line=0,lty=1,lwd=1)
axis(2,at=c(log(0.6),log(1.0), log(2.0), log(3.0)),labels=c("0.6","1.0", "2.0" ,"3.0"), cex.axis=1.05,font.axis=1,las=1,adj=1,line=0,lty=1,lwd=1)
axis(4,at=seq(log(0.6),log(3),length=6),labels=seq(0,(100/(100/30)),length=6),
     cex.axis=1.05,font.axis=1,las=1,hadj=1,line=0,lty=1,lwd=1, mgp=c(5, 2, 0))

mtext(side=1,expression(paste("Sum of Inorganic and Methylated Arsenic, ", symbol("m"), "g/g")),cex=1.10,font=1,line=3.5)

mtext(side=2,"Hazard Ratio",cex=1.10,font=1,las=3,line=2.5)
mtext(side=2,"Percentage of Exposed Participants",cex=1.10,font=1,las=3,line=-25.5)
#text(6.4,0.25,"Percentage of Exposed Participants",cex=1.10,font=1,las=3, srt=270, xpd=NA)

legend(0.3,1.1 ,c("Hazard Ratio", "95% CI"),
       lty=c(1,2), col=c( "black", "black"), cex=0.9,  bty="n")
       
title("Sum of Inorganic and Methylated As")

##################################################################################
#par(mfrow=c(1,3),mar=c(5,5,5,5))

#  	5-B. Plot for inorganic arsenic (Arsenate and Arsenite)

plot(ias,lnHR.ias,type="n",xlim=c(min(ias.breaks),max(ias.breaks)),ylim=c(log(0.6),log(3.00)),xlab="",ylab="",axes=F)    #note that type n is blank

rect(ias.breaks[1:(length(ias.breaks)-1)],rep(log(0.6),length(ias.breaks)-1),				# rect(xleft,ybottom,xright,ytop)
     ias.breaks[2:length(ias.breaks)],log(0.6)+(100/30)*ias.freq*(log(3)-log(0.6)),
     density=-1,border="white",col=grey(0.90),lty=1,lwd=1.5)

lines(ias,lnHR.ias, lty=1,lwd=2, col="black")	
lines(ias,ll.HR.ias,lty=2,lwd=2, col="black")	
lines(ias,ul.HR.ias,lty=2,lwd=2, col="black")			

abline(h=0,lty=1,lwd=1) 

axis(1,at=seq(min(ias.breaks),max(ias.breaks),length=8),labels=round(exp(seq(min(ias.breaks),max(ias.breaks),length=8)),1),
     cex.axis=1.05,font.axis=1,line=0,lty=1,lwd=1)
axis(2,at=c(log(0.6),log(1.0), log(2.0), log(3.0)),labels=c("0.6","1.0", "2.0" ,"3.0"), cex.axis=1.05,font.axis=1,las=1,adj=1,line=0,lty=1,lwd=1)
axis(4,at=seq(log(0.6),log(3),length=6),labels=seq(0,(100/(100/30)),length=6),
     cex.axis=1.05,font.axis=1,las=1,hadj=1,line=0,lty=1,lwd=1, mgp=c(5, 2, 0))

mtext(side=1,expression(paste("Arsenate and Arsenite, ", symbol("m"), "g/g")),cex=1.10,font=1,line=3.5)

mtext(side=2,"Hazard Ratio",cex=1.10,font=1,las=3,line=2.5)
mtext(side=2,"Percentage of Exposed Participants",cex=1.10,font=1,las=3,line=-25.5)
#text(6.4,0.25,"Percentage of Exposed Participants",cex=1.10,font=1,las=3, srt=270, xpd=NA)

legend(0.3,1.1 ,c("Hazard Ratio", "95% CI"),
       lty=c(1,2), col=c( "black", "black"), cex=0.9,  bty="n")
       
title("Arsenate and Arsenite")

#################################################################################

#  	5-C. Plot for MMA

plot(mma,lnHR.mma,type="n",xlim=c(min(mma.breaks),max(mma.breaks)),ylim=c(log(0.6),log(3.00)),xlab="",ylab="",axes=F)    #note that type n is blank

rect(mma.breaks[1:(length(mma.breaks)-1)],rep(log(0.6),length(mma.breaks)-1),  			# rect(xleft,ybottom,xright,ytop)
     mma.breaks[2:length(mma.breaks)],log(0.6)+(100/30)*mma.freq*(log(3)-log(0.6)),
     density=-1,border="white",col=grey(0.90),lty=1,lwd=1.5)

lines(mma,lnHR.mma, lty=1,lwd=2, col="black")	
lines(mma,ll.HR.mma,lty=2,lwd=2, col="black")	
lines(mma,ul.HR.mma,lty=2,lwd=2, col="black")			

abline(h=0,lty=1,lwd=1) 

axis(1,at=seq(min(mma.breaks),max(mma.breaks),length=8),labels=round(exp(seq(min(mma.breaks),max(mma.breaks),length=8)),1),
     cex.axis=1.05,font.axis=1,line=0,lty=1,lwd=1)
axis(2,at=c(log(0.6),log(1.0), log(2.0), log(3.0)),labels=c("0.6","1.0", "2.0" ,"3.0"), cex.axis=1.05,font.axis=1,las=1,adj=1,line=0,lty=1,lwd=1)
axis(4,at=seq(log(0.6),log(3),length=6),labels=seq(0,(100/(100/30)),length=6),
     cex.axis=1.05,font.axis=1,las=1,hadj=1,line=0,lty=1,lwd=1, mgp=c(5, 2, 0))

mtext(side=1,expression(paste("Monomethylarsonate (MMA), ", symbol("m"), "g/g")),cex=1.10,font=1,line=3.5)

mtext(side=2,"Hazard Ratio",cex=1.10,font=1,las=3,line=2.5)
mtext(side=2,"Percentage of Exposed Participants",cex=1.10,font=1,las=3,line=-25.5)
#text(6.4,0.25,"Percentage of Exposed Participants",cex=1.10,font=1,las=3, srt=270, xpd=NA)

#legend(-2.2,1.1 ,c("without iAs adjustment", "with iAs adjustment"),
#      lty=c(1,2), col=c( "black", "red"), cex=0.9,  bty="n")

legend(0.3,1.1 ,c("Hazard Ratio", "95% CI"),
       lty=c(1,2), col=c( "black", "black"), cex=0.9,  bty="n")
       
title("Monomethylarsonate (MMA)")

#################################################################################

#    5-D. Plot for DMA

plot(dma,lnHR.dma,type="n",xlim=c(min(dma.breaks),max(dma.breaks)),ylim=c(log(0.6),log(3.00)),xlab="",ylab="",axes=F)    #note that type n is blank

rect(dma.breaks[1:(length(dma.breaks)-1)],rep(log(0.6),length(dma.breaks)-1),    		# rect(xleft,ybottom,xright,ytop)
     dma.breaks[2:length(dma.breaks)],log(0.6)+(100/30)*dma.freq*(log(3)-log(0.6)),
     density=-1,border="white",col=grey(0.90),lty=1,lwd=1.5)

lines(dma,lnHR.dma, lty=1,lwd=2, col="black")	
lines(dma,ll.HR.dma,lty=2,lwd=2, col="black")	
lines(dma,ul.HR.dma,lty=2,lwd=2, col="black")			

abline(h=0,lty=1,lwd=1) 

axis(1,at=seq(min(dma.breaks),max(dma.breaks),length=8),labels=round(exp(seq(min(dma.breaks),max(dma.breaks),length=8)),1),
     cex.axis=1.05,font.axis=1,line=0,lty=1,lwd=1)
axis(2,at=c(log(0.6),log(1.0), log(2.0), log(3.0)),labels=c("0.6","1.0", "2.0" ,"3.0"), cex.axis=1.05,font.axis=1,las=1,adj=1,line=0,lty=1,lwd=1)
axis(4,at=seq(log(0.6),log(3),length=6),labels=seq(0,(100/(100/30)),length=6),
     cex.axis=1.05,font.axis=1,las=1,hadj=1,line=0,lty=1,lwd=1, mgp=c(5, 2, 0))

mtext(side=1,expression(paste("Dimethylarsinate (DMA), ", symbol("m"), "g/g")),cex=1.10,font=1,line=3.5)

mtext(side=2,"Hazard Ratio",cex=1.10,font=1,las=3,line=2.5)
mtext(side=2,"Percentage of Exposed Participants",cex=1.10,font=1,las=3,line=-25.5)
#text(6.4,0.25,"Percentage of Exposed Participants",cex=1.10,font=1,las=3, srt=270, xpd=NA)

legend(0.3,1.1 ,c("Hazard Ratio", "95% CI"),
       lty=c(1,2), col=c( "black", "black"), cex=0.9,  bty="n")

title("Dimethylarsinate (DMA)")

