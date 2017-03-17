###Plots for OR and HR of reduced eGFR comparing 75th to 25th percentile of arsenic by participant subgroups

##This will include both the OR for cross sectional and the HR for prospective!



#set working directory and read data into R form, then save
setwd("C:/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr")

#at home on mac do this:
setwd("/Users/LauraZheng/Documents/Hopkins Classes Stuff/Strong Heart Study/basic prospective arsenic gfr")

#asplots <-  read.csv("C:/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr/for-plots-in-r.csv", na=" ")
#save(asplots,file="asplots.rda")

asplots <-  read.csv("C:/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr/141211-for-new-forest-plots.csv", na=" ")
save(asplots,file="asplots.rda")


asplots <-  read.csv("/Users/LauraZheng/Documents/Hopkins Classes Stuff/Strong Heart Study/basic prospective arsenic gfr/140829-for-new-forest-plots.csv", na=" ")
save(asplots,file="asplots.rda")



#load data in R form (to bypass the csv crap)
load("C:/Users/lzheng/Documents/strong heart study/basic prospective arsenic gfr/asplots.rda")

load("/Users/LauraZheng/Documents/Hopkins Classes Stuff/Strong Heart Study/basic prospective arsenic gfr/asplots.rda")
# mfrow=c(# rows, # columns)
#par(mfrow=c(1,2),mar=c(4,0.25, 0.25, 0.25))
#par(mfrow=c(1,2),oma=c(2,0,4,0), mar=c(4,3,2,16), lwd=1.10)


# FOR ALL PLOTS  

#Subgroup headings
group <- c(expression(bold("Overall")),
           expression(bold("Sex")),
           expression(bold("Age, years")),
           expression(bold("Study Location")),
           expression(bold("Education")),
           expression(bold(paste(BMI,", ", kg/m^2))),
           expression(bold("Diabetes")),
           expression(bold("Hypertension Treatment")))


#X26egories within the subgroups
label <- c("Men","Women",
           "< 55","55 - 64",expression(paste(symbol("\263")," 65")),  					
           "Arizona","Oklahoma","Dakotas",
           "No HS","Some HS","Completed HS",
           "<30", expression(paste(symbol("\263")," 30")),
           "No", "Yes",
	         "No", "Yes")


# number of people in the subgroups for OR
numOR<-paste(asplots$num_OR)
numOR[c(2,5,9,13,17,20,23)] <- ""    		# change "NA" to blank

#number of people in subgroups for HR
numHR<-paste(asplots$num_HR)
numHR[c(2,5,9,13,17,20,23)] <- ""      	# change "NA" to blank


y1 <- c(23:22,20:18,16:14,12:10,8:7,5:4,2:1)		              # Y-values for OR & HR for each subgroup (ex: < 65 or >65)  (and for N's)
y2 <- c(25,24,21,17,13,9,6,3)								             	      # Y-values for subgroup labels (ex: Age, years)

y3 <- c(25,24,23:22,20:18,16:14,12:10,8:7,5:4,2:1)  	             	# Y-values for OR & HR for each subgroup and overall
y4 <- c(25:1)                                                       # Y-values for OR and HR  (plot blanks too)

#values for OR(95% CI), rounded to 2 digits
#asOR95ci <-paste(format(asplots$as_OR, digits=2, nsmall=2)," (",format(asplots$as_OR_lower,digits=2),"-",format(asplots$as_OR_upper, digits=2),")",sep="")
#asOR95ci[c(2,5,9,13,17,21,25,28,31,34)] <- ""  			# change "NA" to blank

sumOR95ci<-paste(format(asplots$sum_OR, digits=2, nsmall=2)," (",format(asplots$sum_OR_lower,digits=2),"-",format(asplots$sum_OR_upper,digits=2),")",sep="")
sumOR95ci[c(2,5,9,13,17,20,23)] <- ""    		# change "NA" to blank

#values for HR(95% CI), round to 2 digits
#asHR95ci <-paste(format(asplots$as_HR, digits=2, nsmall=2)," (",format(asplots$as_HR_lower,digits=2),"-",format(asplots$as_HR_upper, digits=3),")",sep="")
#asHR95ci[c(2,5,9,13,17,21,25,28,31,34)] <- ""    		# change "NA" to blank

sumHR95ci<-paste(format(asplots$sum_HR, digits=2, nsmall=2)," (",format(asplots$sum_HR_lower,digits=2),"-",format(asplots$sum_HR_upper,digits=3),")",sep="")
sumHR95ci[c(2,5,9,13,17,20,23)] <- ""      	# change "NA" to blank


#p interaction values
#asORp<-paste(format(asplots$total_OR_p, digits=1, nsmall=2))
#asORp[c(1:2,4:5,7:9,11:13,15:17,19:21,23:25,27:28,30:31,33:34, 36:37)] <- ""      	 # change "NA" to blank

sumORp<-paste(format(asplots$sum_OR_p, digits=1, nsmall=2))
sumORp[c(1:2,4:5,7:9,11:13,15:17,19:20,22:23,25)] <- ""                        # change "NA" to blank

#asHRp<-paste(format(asplots$total_HR_p, digits=1, nsmall=2))
#asHRp[c(1:2,4:5,7:9,11:13,15:17,19:21,23:25,27:28,30:31,33:34, 36:37)] <- ""        # change "NA" to blank
#asHRp[c(6)]<-paste("<","0.01",sep="")                                               # substitute "<0.01" for "0.01"

sumHRp<-paste(format(asplots$sum_HR_p, digits=1, nsmall=2))
sumHRp[c(1:2,4:5,7:9,11:13,15:17,19:20,22:23,25)] <- ""                        # change "NA" to blank
#sumHRp[c(6)]<-paste("<","0.01",sep="")


#### Ok we can start plotting!

# mfrow=c(# rows, # columns) this combines the two plots in one window thing

par(mfrow=c(1,2),mar=c(1,1,1,1))
#plots for total arsenic OR
plot(asplots$as_OR,asplots$X26,type='n',log="x",xlab="Odds Ratio (95%CI)", ylab="",axes=F,cex.lab=0.7, xlim=c(0.3, 2))

mtext(side=2, at=40, expression(bold("OR and HR of CKD by total Arsenic, stratified by Participant Characteristics")),las=1,adj=0, cex=1, line=10)  #title

mtext(side=2, at=40,"Characteristics",cex=0.80,font=2,las=1,adj=0,line=16.5)    # top heading for characteristics
mtext(side=2, at=40, "N",font=2, las=1,adj=0, cex=0.80, line=10)              # top heading for N's
mtext(side=2, at=40, "OR (95% CI)",font=2, las=1,adj=0, cex=0.80, line=8)     # for OR(95% CI)
mtext(side=2, at=40, "p-int",font=2, las=1,adj=0, cex=0.80, line=3)           # top heading for P interactions

mtext(side=2,at=y2,group,cex=0.70,font=2,las=1,adj=0,line=16.5)                 # headings for group labels (ie gender, loX26ion)
mtext(side=2,at=y1,label,cex=0.70,font=1,las=1,adj=0,line=15.5)                 # headings for subgroup labels (ie female, arizona)

mtext(side=2,at=y4,numOR,cex=0.70,font=1,las=1,adj=0,line=10)             # gives the N in each subgroup
mtext(side=2,at=y4,asOR95ci,cex=0.70,font=1,las=1,adj=0,line=8)               # gives the OR(95% ci) for each subgroup
mtext(side=2,at=y4,asORp,cex=0.70,font=1,las=1,adj=0,line=3)                  # gives the p-interaction for as*subgroup

#arrows(0.5,y4[33],0.52,y4[33], lty=1, lwd=4, col="red",length=0.01,code=1)   # arrow for UCI that extends beyond x-axis

abline(v=1,lty=1, lwd =2, col = "black")                                        # dashed red line at OR=1
abline(v=c(0.72),lty=3, lwd =2, col = "black")                                 # solid blue line at overall OR
points(asplots$as_OR,asplots$X26,pch=15, col = "black", cex=1.3)                #plots the points
segments(asplots$as_OR_lower,asplots$X26,asplots$as_OR_upper,asplots$X26,lty=1,lwd=2, col = 1)    #plots the lines
axis(1,at=c(0.3,1.0,2),labels=c("0.3","1.0", "2.0"),cex.axis=0.7,font=1,line=1)             #axes


#plots for total arsenic HR
par(mar=c(4,13,4,4))
plot(asplots$as_HR,asplots$X26,type='n',log="x",xlab="Hazard Ratio (95%CI)", ylab="", axes=F,cex.lab=0.7, xlim=c(0.5, 2))

mtext(side=2, at=40, "N",font=2, las=1,adj=0, cex=0.80, line=10)              # top heading for N's
mtext(side=2, at=40, "HR (95% CI)",font=2, las=1,adj=0, cex=0.80, line=8)     # for HR(95% CI)
mtext(side=2, at=40, "p-int",font=2, las=1,adj=0, cex=0.80, line=3)           # top heading for P interactions

mtext(side=2,at=y4,numHR,cex=0.70,font=1,las=1,adj=0,line=10)             # gives the N in each subgroup
mtext(side=2,at=y4,asHR95ci,cex=0.70,font=1,las=1,adj=0,line=8)               # gives the OR(95% ci) for each subgroup
mtext(side=2,at=y4,asHRp,cex=0.70,font=1,las=1,adj=0,line=3)                  # gives the p-interaction for as*subgroup


abline(v=1,lty=1, lwd =2, col = "black")                                        # dashed red line at HR=1
abline(v=c(1.15),lty=3, lwd =2, col = "black")                                 # solid blue line at overall HR
points(asplots$as_HR,asplots$X26,pch=15, col = "black", cex=1.3)              #plots the point
segments(asplots$as_HR_lower,asplots$X26,asplots$as_HR_upper,asplots$X26,lty=1,lwd=2, col = 1)    #plots the lines
axis(1,at=c(0.5,1.0,2.0),labels=c("0.5","1.0","2.0"),cex.axis=0.7,font=1,line=1)             #axes







#############################################################################################################


#### Plot for sum of inorganic and methylated arsenic

# mfrow=c(# rows, # columns) this combines the two plots in one window thing

par(mfrow=c(1,2),mar=c(4,17,4,0))
#plots for total arsenic OR
plot(asplots$sum_OR,asplots$X26,type='n',log="x",xlab="Odds Ratio (95%CI)", ylab="",axes=F,cex.lab=0.7, xlim=c(0.4, 1.5))

mtext(side=2, at=27.5, expression(bold("OR and HR of CKD by sum of inorganic and methylated arsenic, stratified by Participant Characteristics")),las=1,adj=0, cex=1, line=16)  #title

mtext(side=2, at=26,"Characteristics",cex=0.80,font=2,las=1,adj=0,line=16.5)     # top heading for characteristics
mtext(side=2, at=26, "N",font=2, las=1,adj=0, cex=0.80, line=10)                 # top heading for N's
mtext(side=2, at=26, "OR (95% CI)",font=2, las=1,adj=0, cex=0.80, line=8)        # for OR(95% CI)
mtext(side=2, at=26, "p-int",font=2, las=1,adj=0, cex=0.80, line=3)              # top heading for P interactions

mtext(side=2,at=y2,group,cex=0.70,font=2,las=1,adj=0,line=16.5)                  # headings for group labels (ie gender, loX26ion)
mtext(side=2,at=y1,label,cex=0.70,font=1,las=1,adj=0,line=15.5)                  # headings for subgroup labels (ie female, arizona)

mtext(side=2,at=y4,numOR,cex=0.70,font=1,las=1,adj=0,line=10)                    # gives the N in each subgroup
mtext(side=2,at=y4,sumOR95ci,cex=0.70,font=1,las=1,adj=0,line=8)                 # gives the OR(95% ci) for each subgroup
mtext(side=2,at=y4,sumORp,cex=0.70,font=1,las=1,adj=0,line=3)                    # gives the p-interaction for as*subgroup

#arrows(0.5,y4[33],0.52,y4[33], lty=1, lwd=4, col="red",length=0.01,code=1)      # arrow for UCI that extends beyond x-axis

abline(v=1,lty=1, lwd =2, col = "black")                                         # dashed red line at OR=1
abline(v=c(0.68),lty=3, lwd =2, col = "black")                                   # solid blue line at overall OR
points(asplots$sum_OR,asplots$X26,pch=15, col = "black", cex=1.3)                #plots the points
segments(asplots$sum_OR_lower,asplots$X26,asplots$sum_OR_upper,asplots$X26,lty=1,lwd=2, col = 1)    #plots the lines
axis(1,at=c(0.4,1.0,1.5),labels=c("0.4","1.0", "1.5"),cex.axis=0.7,font=1,line=1)             #axes

#plots forHR
par(mar=c(4,13,4,4))
plot(asplots$sum_HR,asplots$X26,type='n',log="x",xlab="Hazard Ratio (95%CI)", ylab="", axes=F,cex.lab=0.7, xlim=c(0.8, 1.7))

mtext(side=2, at=26, "N",font=2, las=1,adj=0, cex=0.80, line=11)                 # top heading for N's
mtext(side=2, at=26, "HR (95% CI)",font=2, las=1,adj=0, cex=0.80, line=9)        # for HR(95% CI)
mtext(side=2, at=26, "p-int",font=2, las=1,adj=0, cex=0.80, line=4)              # top heading for P interactions

mtext(side=2,at=y4,numHR,cex=0.70,font=1,las=1,adj=0,line=11)                # gives the N in each subgroup
mtext(side=2,at=y4,sumHR95ci,cex=0.70,font=1,las=1,adj=0,line=9)                 # gives the OR(95% ci) for each subgroup
mtext(side=2,at=y4,sumHRp,cex=0.70,font=1,las=1,adj=0,line=4)                    # gives the p-interaction for as*subgroup

arrows(1.70,y4[6],1.75,y4[6], lty=1, lwd=2, col="black",length=0.05,code=2)      # arrow for UCI that extends beyond x-axis
arrows(1.70,y4[14],1.75,y4[14], lty=1, lwd=2, col="black",length=0.05,code=2)      # arrow for UCI that extends beyond x-axis


abline(v=1,lty=1, lwd =2, col = "black")                                         # dashed red line at HR=1
abline(v=c(1.20),lty=3, lwd =2, col = "black")                                   # solid blue line at overall HR
points(asplots$sum_HR,asplots$X26,pch=15, col = "black", cex=1.3)                 #plots the point
segments(asplots$sum_HR_lower,asplots$X26,asplots$sum_HR_upper,asplots$X26,lty=1,lwd=2, col = 1)    #plots the lines
axis(1,at=c(0.8,1.0,1.7),labels=c("0.8","1.0","1.7"),cex.axis=0.7,font=1,line=1)             #axes

