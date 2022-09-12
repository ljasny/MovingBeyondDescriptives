##1 Getting Started
##Statnet has recently been updated and this code is based on the current latest version of R (4.1.0) and Statnet (4.0.0)
library(statnet)

##The sna and network packages are part of statnet, but might need to be loaded separately to attach data
#library(network)
#library(sna)

###########################################################
###########################################################
##2 Node Level Indices, Network Covariates, and Network Regression
##2.1 Network visualization and basic descriptives

data(emon) # Load Drabek et al. data

# Begin by examining the 'emon' dataset
?emon
class(emon)
class(emon[[1]])
emon[[1]]
emon[[1]]%v%"vertex.names" # Display vertex names

# Extract ties from the Cheyenne EMON communicating at least "every few hours"
g<-as.sociomatrix(emon[[1]],"Frequency") # Need to get the frequency info
g<-symmetrize((g>0)&(g<4)) # Note the reverse coding! The frequency of communication is weighted data

#Get some potential covariates
drs<-emon[[1]]%v%"Decision.Rank.Score" # Get decision rank (see man page)
crs<-emon[[1]]%v%"Command.Rank.Score" # Get command rank

# Plot Cheyenne EMON communications
gplot(emon[[1]])
gplot(emon[[1]], label.cex=0.5, label.col=4, label=network.vertex.names(emon[[1]])) # Basic display, with labels

#Calculate some basic centrality measures
deg<-degree(g,gmode="graph")
bet<-betweenness(g,gmode="graph")
clo<-closeness(g,gmode="graph")

#Raw correlations
cor(cbind(deg,bet,clo),cbind(drs,crs))

#Classical tests (using asymptotic t distribution)
cor.test(deg,drs)
cor.test(bet,drs)
cor.test(clo,drs)

###########################################################
#2.2 Testing correlations
#Lets build a permutation test
deg
drs
c.obs<-cor(deg,drs)
c.obs

#Permute one of the data sets
sample(drs)
cor(deg,sample(drs))

#write a for loop to permute one of the data sets many times
c.rep<-vector(length=1000)

for(n in 1:1000){
  c.rep[n]<-cor(deg,sample(drs))
}

#look at what we've created
c.rep
hist(c.rep)

#compare to empirical data
abline(v=c.obs,col="red")

#put it all into a function!
perm.cor.test<-function(x,y,n=1000){ #Define a simple test function
  c.obs<-cor(x,y)
  c.rep<-vector(length=n)
  for(i in 1:n)
    c.rep[i]<-cor(x,sample(y))
  c.rep
}

#look at the results
new.c.rep<-perm.cor.test(deg,drs)
hist(new.c.rep)
abline(v=c.obs,col="red")

#calculate quantiles
mean(new.c.rep>=c.obs)
mean(new.c.rep<=c.obs)
mean(abs(new.c.rep)>=abs(c.obs))

#make the output display these quantiles
perm.cor.test<-function(x,y,niter=1000){ #Define a simple test function
  c.obs<-cor(x,y)
  c.rep<-vector()
  for(i in 1:niter)
    c.rep[i]<-cor(x,sample(y),use="complete.obs")
  cat("Vector Permutation Test:\n\tObserved correlation: ",c.obs,"\tReplicate quantiles
      (number of iterations = ",n,")\n",sep="")
  cat("\t\tPr(rho>=obs):",mean(c.rep>=c.obs),"\n")
  cat("\t\tPr(rho<=obs):",mean(c.rep<=c.obs),"\n")
  cat("\t\tPr(|rho|>=|obs|):",mean(abs(c.rep)>=abs(c.obs)),"\n")
  invisible(list(obs=c.obs,rep=c.rep))
}

#examine the fancy output
fancy<-perm.cor.test(deg,drs)
attributes(fancy)
fancy$obs
fancy$rep

#For more information....
?cor.test
?t.test
?sample

###########################################################
#2.3 Using NLIs as regression covariates

pstaff<-emon[[1]]%v%"Paid.Staff" # Get more EMON covariates
vstaff<-emon[[1]]%v%"Volunteer.Staff"
govt<-((emon[[1]]%v%"Sponsorship")!="Private")

#Simple model: decision rank is linear in size, degree, and govt status
mod<-lm(drs~deg+pstaff+vstaff+govt)
summary(mod)
anova(mod) #Some useful lm tools
AIC(mod)

#Try with alternative measures....
mod2<-lm(drs~bet+pstaff+vstaff+govt) #Betweenness
summary(mod2)
mod3<-lm(drs~clo+pstaff+vstaff+govt) #Closeness
summary(mod3)
AIC(mod,mod2,mod3) #Closeness wins!

#For more information....
?lm
?anova
?AIC

###########################################################
#2.4 Graph correlation and bivariate QAP
# Remember the Florentine families data ## good to use also with attribute data
data(florentine)
gplot(flobusiness) # Examine business ties
gplot(flomarriage) # Examine marriage ties

# Could the two be related?
par(mfrow=c(1,2))
coords<-plot(flobusiness,label=flobusiness%v%"vertex.names",label.cex=.5,pad=1)
title("Business Ties")
plot(flomarriage,label=flomarriage%v%"vertex.names",label.cex=.5,pad=1,coord=coords)
title("Marriage Ties")
par(mfrow=c(1,1))

# Let's try a graph correlation
gcor(flobusiness,flomarriage)

# Why can't we use our previous permutation test?
# instead, use rmperm
# take a look
rmperm

par(mfrow=c(1,2))
gplot(flobusiness[,],label=flobusiness%v%"vertex.names",label.cex=.5,pad=1,coord=coords)
title("Business Ties")
gplot(rmperm(flobusiness),label=flobusiness%v%"vertex.names",label.cex=.5,pad=1,coord=coords)
title("Permuted Business Ties")
par(mfrow=c(1,1))

#now look at qaptest
qaptest

#and use it
qt<-qaptest(list(flobusiness,flomarriage),gcor,g1=1,g2=2)
summary(qt) # Examine the results
plot(qt) # Plot the QAP distribution

# Testing against covariate effects
wealth<-sapply(flomarriage%v%"wealth",rep,network.size(flomarriage))
wealth
wealthdiff<-abs(outer(flomarriage%v%"wealth",flomarriage%v%"wealth","-"))
wealthdiff
qt1<-qaptest(list(flomarriage,wealth),gcor,g1=1,g2=2)
qt2<-qaptest(list(flomarriage,wealthdiff),gcor,g1=1,g2=2)
summary(qt1) # Do wealthy families have more ties?
summary(qt2) # Is there a wealth difference effect?

# For more information....
?qaptest
?gcor
?outer
?sapply
?rep

###########################################################
#2.5 Network Regression

# Combine the previous tests (takes a while to perform QAP test)
marriagefit<-netlm(flomarriage,list(flobusiness,wealth,wealthdiff))
summary(marriagefit) # Examine the results


# Another example: we begin by preparing the response variable. We will use the Cheyenne
# EMON in valued form, but need to recode the frequency data
Y<-as.sociomatrix(emon[[1]], "Frequency") # Extract frequencies
Y[Y>0]<-5-Y[Y>0] # Now, higher -> more frequent

# Extract some vertex attributes
crk<-emon[[1]]%v% "Command.Rank.Score" # Command rank
spon<-emon[[1]]%v%"Sponsorship" # Type of organization

# Form some predictor matrices (X variables)
Xcr<-sapply(crk,rep,length(crk)) # Column effects for command rank
Xcr
Xsp<-outer(spon,spon,"!=") # Dyadic effect for type difference
Xsp

# Fit the model (takes a while to perform QAP test)
cmfit<-netlm(Y,list(Xcr,Xsp))
summary(cmfit) # Examine the results


cmfitB<-netlm(Y,list(Xcr,Xsp),nullhyp="classical")
summary(cmfitB) # In this example, pretty similar

# For more information....
?outer
?netlm

#######################################
##2.6 Network Autocorrelation models

#Network Autocorrelation Models

#We need another library
#install.packages("numDeriv")
library(numDeriv)

##and to load some data. 
##Make sure your working directory is set correctly
load("movingBeyondDescriptivesData.Rdata")

rownames(parishNet)
which(parishNet["Rapides",]==1)

#Examine the OLS regression predicting percentage vote for Obama by percentage Black residents and average weekly wage

modLM<-lm(parishAttr[,1]~parishAttr[,2]+parishAttr[,3])
summary(modLM)

#create some simple weight matrices
wMat<-t(sapply(1:nrow(parishNet),function(x){parishNet[x,]/sum(parishNet[x,])}))
wRandMat<-rgraph(64)

#Run the different models with the empirical weights and look at the results 

modAR<-lnam(parishAttr[,1],x=parishAttr[,2:3],W1=wMat) ##put in as a matrix what your covariates are
summary(modAR)

modMA<-lnam(parishAttr[,1],x=parishAttr[,2:3],W2=wMat)
summary(modMA)

modARMA<-lnam(parishAttr[,1],x=parishAttr[,2:3],W1=wMat,W2=wMat)
summary(modARMA)

modRandAR<-lnam(parishAttr[,1],x=parishAttr[,2:3],W1=wRandMat)
summary(modRandAR)

modWIntercept<-lnam(parishAttr[,1],x=cbind(1,parishAttr[,2:3]),W1=wMat)
summary(modWIntercept)


###########################################################
###########################################################
##3 Simple univariate conditional uniform graph tests
# The cug.test function provides a simple interface for univariate CUG tests.
# Let's try testing some data on trade in complex manufactured goods to see if overall
# activity (density) is greater then would be expected from size alone.
cug.test(ctrade,gden) # Call cug.test with the gden (density) function

# Is there more reciprocity than density would suggest? Let's see.
cug.test(ctrade,grecip,cmode="edges") # Conditioning on edges, calling grecip

# Given biases in density and reciprocity, we might look to see if there is a
# bias in transitivity, too. This time, let's condition on all of the above.
cug.test(ctrade,gtrans,cmode="dyad") # Conditioning on dyad census

# We are not limited to simple commands. Let's try indegree centralization:
ct<-cug.test(ctrade,centralization,cmode="dyad",FUN.arg=list(FUN=degree, cmode="indegree")) 
# Note that we here must pass not only arguments to
# cug.test, but also to centralization and degree!

ct # Print the result
plot(ct) # Can also plot it!

###############################
##3.2 Using CUG to produce reference distributions

mod1<-netlm(ctrade,list(trade[1,,],trade[3,,],trade[4,,],trade[5,,]),nullhyp = "classical")

mod2<-netlm(ctrade,list(trade[1,,],trade[3,,],trade[4,,],trade[5,,]),nullhyp = "qap")

mod3<-netlm(ctrade,list(trade[1,,],trade[3,,],trade[4,,],trade[5,,]),nullhyp = "cugden")

