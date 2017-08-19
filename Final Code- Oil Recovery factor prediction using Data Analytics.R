#Thesis Project
#Name: Gowtham Talluru

#_____________________________________________________________________________________________________________________
##############################################################################################################################
#required packages for the Project
library(robustbase)                     # Adjusted Boxplot  
library(EnvStats)                       # box-cox transformation
library(rpart)                          # for decision tree modeling
library(party)                          # for visualizing trees
library(partykit)                       # for visualizing trees
library(rattle)      		                # fancy tree plot
library(rpart.plot)			                # enhanced tree plots
library(caret)
library(RANN)
library(pls)
library(corrplot)
library(car)           
library(DescTools)
library(psych)
library(Boruta)
library(XLConnect)
library(randomForest)
library(neuralnet)
library(nnet)
library(glmnet)           #Penalized regression
library(flexclust)
library(readxl) 
library(class)            #For k nearest neighbours
library(MASS)             #Robust Regression
library(Hmisc)            #For dataframe histogram
library(e1071)            #For measuring skewness
#_____________________________________________________________________________________________________________________
##############################################################################################################################

# Sequence of Analysis
# 1. Initial Data Analysis (Loading, Str, Dim,complete case, missing if any)
# 2. Formatting the predictors for different representations and converting them to numeric form
# 3. Dividing data for oil and gas reservoir (Based on Sand type and GOR/ removing water sands)
# 4. Releveling of categorical variables
# 5. Outlier Identification 
# 6. Splitting and Plotting of numeric and categorical variable 
# 7. Normalization of skewed numerical data (If it is required on both test/train data)
# 8. Subselecting predictors for each reservoir type (Both test/train data)
# 9. Variable Importance determination on train data: Decision tree/ Boruta/ Anova
# 10.Dividing data between test and train for each reservoir type
# 11.Model fitting for prediction
# 12.Cross Validation





#=====================================================================================================================
#---------------------------------------------------------------------------------------------------------------------
#                             ************** DATA UNDERSTANDING ****************
#---------------------------------------------------------------------------------------------------------------------
#=====================================================================================================================

#####################################################################################
#                         Initial Evaluation of the data set                       #
#####################################################################################
#------------------------------------------------------------------------------------

# Loading Data
setwd("C:\\Users\\Gowtham\\Desktop\\R-2")
rawdata <- read.csv("rawdata.csv")
# Checking the Dimensions of raw Data
dim(rawdata)                      
# 13289 x 82

# Checking Structure of Data
str (rawdata)

# 1st Filtering: Complete Cases
# Checking No. of complete cases available and checkig its dimensions
rawDataComlCases<-rawdata[complete.cases(rawdata),]
dim(rawDataComlCases)
# 13288 x 82

#=====================================================================================================================
#---------------------------------------------------------------------------------------------------------------------
#                                 ************** DATA PREPARATION ****************
#---------------------------------------------------------------------------------------------------------------------
#=====================================================================================================================


#####################################################################################
#                       Changing predictors to right format                         #
#####################################################################################
#------------------------------------------------------------------------------------

rawDataComlCases$P_CUMCOIL<-as.numeric(gsub(",","",rawDataComlCases$P_CUMCOIL))  # Np
rawDataComlCases$P_RECOIL<-as.numeric(gsub(",","",rawDataComlCases$P_RECOIL))    # Proved recoverable oil
rawDataComlCases$ORECO<-as.numeric(gsub(",","",rawDataComlCases$ORECO))          # Oil Reservoirs proved recoverable oil
rawDataComlCases$P_RECBOE<-as.numeric(gsub(",","",rawDataComlCases$P_RECBOE))    # Proved Recoverable BOE
rawDataComlCases$P_CUMBOE<-as.numeric(gsub(",","",rawDataComlCases$P_CUMBOE))    # Cumulative BOE produced
rawDataComlCases$DISCBOE<-as.numeric(gsub(",","",rawDataComlCases$DISCBOE))      # Discoverable BOE
rawDataComlCases$SS<-as.numeric(gsub(",","",rawDataComlCases$SS))                # Subsea depth
rawDataComlCases$TAREA<-as.numeric(gsub(",","",rawDataComlCases$TAREA))          # Total Area
rawDataComlCases$OAREA<-as.numeric(gsub(",","",rawDataComlCases$OAREA))          # Oil Area
rawDataComlCases$PERMEABILITY<-as.numeric(gsub(",","",rawDataComlCases$PERMEABILITY))      # Permeability
rawDataComlCases$PI<-as.numeric(gsub(",","",rawDataComlCases$PI))                # Initial Reservoir Pressure
rawDataComlCases$WDEP<-as.numeric(gsub(",","",rawDataComlCases$WDEP))            # Water Depth

#####################################################################################
# Separating Oil and Gas Reservoirs and applying filters 
#####################################################################################

# 2nd Filter:  Based on Sand Type
#------------------------------------------------------------------------------------
 
oilRes<-rawDataComlCases[rawDataComlCases$SD_TYPE=="O" | rawDataComlCases$SD_TYPE==0
                         |rawDataComlCases$SD_TYPE=="B",]


# 3rd Filter: Based on Gas Oil ratio, generally GOR < 30 Mcf/bbl is Oil Reservoir
#------------------------------------------------------------------------------------
oilRes<-oilRes[oilRes$GOR<30,]
dim(oilRes)
# 4313 x 82

# 4th Filter:Dropping Reservoirs in which very few categories are present
#------------------------------------------------------------------------------------
# Checking distribution in Play type
table(oilRes$PLAY_TYPE)
# Subsetting play types with maximum observations
oilRes<-oilRes[(oilRes$PLAY_TYPE!='B1')&(oilRes$PLAY_TYPE!='R1')&    
                 (oilRes$PLAY_TYPE!='X2'),]
dim(oilRes)
# 4244 x 82

#--------------#
# Checking distribution in Formation structure
table(oilRes$FSTRU)
# Subsetting Formation structure with maximum observations
oilRes<-oilRes[(oilRes$FSTRU!='F')&(oilRes$FSTRU!='G')&
                 (oilRes$FSTRU!='H')&(oilRes$FSTRU!='I'),]
dim(oilRes)
# 4011 x 82

#--------------#
# Checking distribution in Drive 
table(oilRes$DRIVE)
# Subsetting Drive with maximum observations
oilRes<-oilRes[(oilRes$DRIVE!='0')&(oilRes$DRIVE!='GCP')&
                 (oilRes$DRIVE!='SLG')&(oilRes$DRIVE!='UNK'),]
dim(oilRes)
# 3724 x 82

#--------------#
# Checking distribution in Reservoir type
table(oilRes$RES_TYPE)  
#    N    S    U                 # Even Distribution
#   701  721 2302 

#--------------#
# Checking the distributions in Refined Data set
boxplot(ORF~DRIVE, data=oilRes)
boxplot(ORF~FSTRU, data=oilRes)
boxplot(ORF~PLAY_TYPE, data=oilRes)
boxplot(ORF~RES_TYPE, data=oilRes)

dim(oilRes)
# 3724 x 82

# 5th Filter: Dropping reservoirs having zero OIP,ORF and zero BHCOMP
#------------------------------------------------------------------------------------
oilRes<-oilRes[oilRes$OIP!=0,]
oilRes<-oilRes[oilRes$BHCOMP!=0,]
oilRes<-oilRes[oilRes$PERMEABILITY!=0,]
oilRes<-oilRes[oilRes$ORF!=0,]
oilRes$GOR<-ifelse(oilRes$GOR==0, 1,oilRes$GOR)
dim(oilRes)
# 3342 x  82

# 6th Filter: Taking Reservoirs with completed production
#------------------------------------------------------------------------------------
oilRes<-oilRes[oilRes$P_CUMCOIL>(0.8*oilRes$P_RECOIL),]
dim(oilRes)
#3038 x 82

# 7th Filter: Discrepancy between reported Oil Recovery Factor (ORF) and Calculated Recovery
#------------------------------------------------------------------------------------
# Cross checking ORF and Np/OIP
rf<-oilRes$P_CUMCOIL/oilRes$OIP
plot(rf, oilRes$ORF, xlim=c(0,1), ylim=c(0,1),
     main="Recovery Factor(RF) reporter vs Calculated", col="red",
     xlab="RF Calculated", ylab="RF reported")

#--------------#
# There is a lot of discripency between two values
# Subselecting the reservoirs whose ORF is near to Np/OIP
oilResClean<-oilRes[(abs(rf-oilRes$ORF)<0.05),]

dim(oilResClean)
# 1462 x 82

#--------------#
# Double check by Ploting the same graph for cleaned oil Reservoirs
rClean<-oilResClean$P_CUMCOIL/oilResClean$OIP
plot(rClean, oilResClean$ORF, xlim=c(0,1), ylim=c(0,1),
     main="Recovery Factor(RF) reporter vs Calculated", col="red",
     xlab="RF Calculated", ylab="RF reported")

dim(oilResClean)
# 12517 x 82
write.csv(oilRes, " Subset.csv")
################################################################################
################################################################################
################################################################################
#*******************************************************************************
#*******************************************************************************
#Taking in the dimless parameters calculated in Excel
oilResClean <- read_excel("Dimless_Calculation.xlsm")

oilResClean$CHRONOZONE2<-oilResClean$CHRONOZONE

#1462 x 94
#####################################################################################
#                        Releveling Categorical Variables                          #
#####################################################################################

# Reducing the number of factors in CHRONOZONE
#------------------------------------------------------------------------------------
oilResClean$CHRONOZONE2[(oilResClean$CHRONOZONE=="MML")|
                          (oilResClean$CHRONOZONE=="MUL")|
                          (oilResClean$CHRONOZONE=="MLM")|
                          (oilResClean$CHRONOZONE=="mMMM")|
                          (oilResClean$CHRONOZONE=="MUM")|
                          (oilResClean$CHRONOZONE=="MMM")|
                          (oilResClean$CHRONOZONE=="MLL")|
                          (oilResClean$CHRONOZONE=="MUU")|
                          (oilResClean$CHRONOZONE=="KUL")|
                          (oilResClean$CHRONOZONE=="MLU")]="MIOCENE"

#oilResClean$CHRONOZONE2[(oilResClean$CHRONOZONE=="PL")]="PLIO_LWR"

oilResClean$CHRONOZONE2[(oilResClean$CHRONOZONE=="PU")|
                          (oilResClean$CHRONOZONE=="PU-PL")]="PLIO_UPR"


oilResClean$CHRONOZONE2[(oilResClean$CHRONOZONE=="PLU-LL")|
                          (oilResClean$CHRONOZONE=="PLM")|
                          (oilResClean$CHRONOZONE=="PL")|
                          (oilResClean$CHRONOZONE=="PLU")|
                          (oilResClean$CHRONOZONE=="PLL")]="PLEISTOCENE"

oilResClean$CHRONOZONE2<-as.factor(oilResClean$CHRONOZONE2)


# Merging reservoirs with index 0 to UNK
#oilResClean$DRIVE[oilResClean$DRIVE=='0']<-'UNK'


# Merging Gas cap drive reservoirs to Combination drive
# GCP ----> COM
oilResClean$DRIVE[oilResClean$DRIVE=='GCP']<-'COM'
table(oilResClean$DRIVE)



#=====================================================================================================================
#---------------------------------------------------------------------------------------------------------------------
#                       **************  DATA UNDERSTANDING CONTINUED ****************
#---------------------------------------------------------------------------------------------------------------------
#=====================================================================================================================




#####################################################################################
#                Variable Selection                 #
#####################################################################################

# There are a lot of duplicate variables in the data
# For example there are 4 different identifiers for each hydrocarbon reservoir
# The following parameters are selected for further model building

# Final selecton of predictors
#------------------------------------------------------------------------------------
#[,c("FCLASS","FSTRU","PLAY_TYPE","FTRAP1","DRIVE","RES_TYPE")]
rCleanSelect<-oilResClean[,c("SAND_NAME",
                             "FSTRU","PLAY_TYPE",
                             "CHRONOZONE2","SPGR",
                             "DRIVE","RES_TYPE","POROSITY", "SW",
                             "BHCOMP","Mobility_Ratio_endpoint",
                             "Nalpha","RL","Np","Ng","Npc","ORF","Dev_factor","Hetro_factor")]

#rCleanSelect$FCLASS<-as.factor(rCleanSelect$FCLASS)
rCleanSelect$FSTRU<-as.factor(rCleanSelect$FSTRU)
rCleanSelect$PLAY_TYPE<-as.factor(rCleanSelect$PLAY_TYPE)
rCleanSelect$CHRONOZONE2<-as.factor(rCleanSelect$CHRONOZONE2)
rCleanSelect$DRIVE<-as.factor(rCleanSelect$DRIVE)
rCleanSelect$RES_TYPE<-as.factor(rCleanSelect$RES_TYPE)
rCleanSelect$SAND_NAME<-as.factor(rCleanSelect$SAND_NAME)
#rCleanSelect$FTRAP1<-as.factor(rCleanSelect$FTRAP1)


#  Naming rows with sand names
rownames(rCleanSelect)<-rCleanSelect$SAND_NAME
rCleanSelect$SAND_NAME<-NULL



# Step 2: Using Correlation Matrix for selected feature engineering and variable selection
#------------------------------------------------------------------------------------
oilResClean_int <- rCleanSelect[,sapply(rCleanSelect, is.numeric)]      # Separating the integer variables

c <-cor(oilResClean_int)

# Correlation Matrix
corrplot.mixed(c,lower = "number", upper = "circle",order="AOE", diag="n", tl.pos = "lt", 
               tl.cex = 0.7,cl.cex = 0.7, number.cex = 0.6, na.rm=TRUE )

oilResClean_int<-NULL                                   #Removing extra variable
c<-NULL                                                 #Removing extra variable

#-------------------------------------------------------------------------------
#Removing absurd/erroneous values from data set
#-------------------------------------------------------------------------------
summary(rCleanSelect$Mobility_Ratio_endpoint)
hist(rCleanSelect$Mobility_Ratio_endpoint, col="blue", 
     main ="End Point Mobility Ratio", xlab="")


summary(rCleanSelect$Nalpha)
hist(rCleanSelect$Mobility_Ratio_endpoint, col="red", 
     main ="Histogram N-alpha", xlab="N-alpha")


summary(rCleanSelect$Np)
hist(rCleanSelect$Mobility_Ratio_endpoint, col="green", 
     main ="Histogram Np", xlab="Np")



#Removing reservoirs having absurd values
rCleanSelect<-rCleanSelect[rCleanSelect$Mobility_Ratio_endpoint<10, ]
rCleanSelect<-rCleanSelect[rCleanSelect$Nalpha<5,]
rCleanSelect<-rCleanSelect[rCleanSelect$Np<15,]
rCleanSelect<-rCleanSelect[rCleanSelect$Hetro_factor<10,]


#####################################################################################
#   Seperating Numeric and Categorical predictors in the Cleaned Data Set           #
#####################################################################################

l<-sapply(rCleanSelect,function(x)!is.factor(x))          #Numeric variables in rCleanSelect
numericVariables<-names(rCleanSelect)[which(l=="TRUE")]   #Numeric variables in rCleanSelect
factorVariables<-names(rCleanSelect)[which(l=="FALSE")]   #Factor variables in rCleanSelect 


GoM_factors<-rCleanSelect[factorVariables]       #Seperating Categorical variables
GoM_factors$ORF<-rCleanSelect$ORF                #Adding outcome varibale to factors
dim(GoM_factors)
# 3038   7
GoM_numeric<-rCleanSelect[numericVariables]      #Seperating Numeric variables
dim(GoM_numeric)
data.frame(GoM_numeric2)
# 3038   15

################################################################################
#Inspecting each predictor
################################################################################

# 1.Summary of Specific gravity
summary(rCleanSelect$SPGR)
hist(rCleanSelect$SPGR)
GoM_numeric2$SPGR<-((rCleanSelect$SPGR-mean(rCleanSelect$SPGR))/sd(rCleanSelect$SPGR))
hist(GoM_numeric2$SPGR)


# 2.Summary of Porsity
summary(rCleanSelect$POROSITY)
hist(rCleanSelect$POROSITY)
GoM_numeric2$POROSITY<-((rCleanSelect$POROSITY-mean(rCleanSelect$POROSITY))/sd(rCleanSelect$POROSITY))
hist(GoM_numeric2$POROSITY)

# 3. Summary of Sw
summary(rCleanSelect$SW)
hist(rCleanSelect$SW)
GoM_numeric2$SW<-((rCleanSelect$SW-mean(rCleanSelect$SW))/sd(rCleanSelect$SW))
hist(GoM_numeric2$SW)


# 4. Summary of BHCOMP (Number of completions)
summary(rCleanSelect$BHCOMP)
hist(rCleanSelect$BHCOMP)
temp<-log(rCleanSelect$BHCOMP)
GoM_numeric2$BHCOMP<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$BHCOMP)

 
# 5. Summary of Mobility_Ratio_endpoint
temp<-rCleanSelect$Mobility_Ratio_endpoint
summary(temp)
hist(temp)
GoM_numeric2$Mobility_Ratio_endpoint<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Mobility_Ratio_endpoint)


# 6. Summary of Nalpha
temp<-rCleanSelect$Nalpha
summary(temp)
hist(temp)
hist(log(temp))
temp<-log(temp)
GoM_numeric2$Nalpha<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Nalpha)


# 7. Summary of RL
temp<-rCleanSelect$RL
summary(temp)
hist(temp)
hist(log(temp))
temp<-log(temp)
GoM_numeric2$RL<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$RL)

# 7. Summary of Np
temp<-rCleanSelect$Np
summary(temp)
hist(temp)
GoM_numeric2$Np<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Np)

#8. Summary of Ng
temp<-rCleanSelect$Ng
summary(temp)
hist(temp)
hist(log(temp))
temp<-log(temp)
GoM_numeric2$Ng<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Ng)

#9. Summary of Npc
temp<-rCleanSelect$Npc
summary(temp)
hist(temp)
hist(log(temp))
temp<-log(temp)
GoM_numeric2$Npc<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Npc)

#10. Summary of Dev_factor
temp<-rCleanSelect$Dev_factor
summary(temp)
hist(temp)
hist(log(temp))
temp<-log(temp)
GoM_numeric2$Dev_factor<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Dev_factor)

#11. Summary of Hetro_factor
temp<-rCleanSelect$Hetro_factor
summary(temp)
hist(temp)
GoM_numeric2$Hetro_factor<-(temp-mean(temp))/sd(temp)
hist(GoM_numeric2$Hetro_factor)


GoM_numeric2$ORF<-rCleanSelect$ORF
GoM_numeric2<-as.data.frame(GoM_numeric2)

# For all the models (except DT and RF), it behooves to scale numeric attributes 
# and introduce dummy varibles. Therefore, following steps were performed. 

################################################################################
#  Pre-Processing Numeric (Scaling) and Factor Variables (Dummy Variables)     #
################################################################################


# Step 2: Introducing dummy variables for factor variables
#------------------------------------------------------------------------------------
# Transforming dummy variables
GoM_factors$ORF<-NULL

dm<-dummyVars("~.",data=GoM_factors)
GoM_factors_Dummy<-data.frame(predict(dm, newdata=GoM_factors))
dim(GoM_factors_Dummy)
#2524 * 22


#################################################################################
#     Creating various data sets (Original, processed, processed w/o outliers)
#################################################################################
GoM_original<-rCleanSelect

GoM_processed<-cbind(GoM_factors_Dummy, GoM_numeric2)

#Dataset without outliers
GoM_proc_nooutliers<-GoM_processed
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_processed$SPGR)<=3,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$POROSITY)<=3,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$SW)<=3.5,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Mobility_Ratio_endpoint)<=3.5,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Nalpha)<=3,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$RL)<=4,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Np)<=3.5,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Ng)<=4,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Npc)<=4,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Dev_factor)<=4,]
GoM_proc_nooutliers<-GoM_proc_nooutliers[abs(GoM_proc_nooutliers$Hetro_factor)<=4,]

dim(GoM_proc_nooutliers)

#####################################################################################
#                  Splitting into Test and Train Data Set                           #
#####################################################################################
#Splitting original data and processed data

# Setting seed for reproducibility
set.seed(1)

# Data is split by keeping same distribution of ORF in train and test data
train<-createDataPartition(GoM_original$ORF, p=0.80,list=FALSE)

GoM_processed$SPGR<-NULL
GoM_processed$RL<-NULL

GoM_proc_nooutliers$SPGR<-NULL
GoM_proc_nooutliers$RL<-NULL
#Train and Test dataset
GoMTrain_original<-GoM_original[train,]
GoMTrain_processed<-GoM_processed[train,]


GoMTest_original<-GoM_original[-train,]
GoMTest_processed<-GoM_processed[-train,]

#Splitting dataset with out outliers
set.seed(1)
train2<-createDataPartition(GoM_proc_nooutliers$ORF, p=0.80,list=FALSE)

GoMTrain_proc_nooutliers<-GoM_proc_nooutliers[train2,]
GoMTest_proc_nooutliers<-GoM_proc_nooutliers[-train2,]


#Plots   for      ppt
#Data distribution
par(mfrow=c(1,3))
hist(GoM_original$POROSITY, col="red", main="GoM_Original", xlab="POROSITY")
hist(GoM_processed$POROSITY, xlim=c(-4,4), main="GoM_Processed", col="red", 
     xlab="POROSITY")
abline(v=c(-3,3), col="blue",lwd=2)
hist(GoM_proc_nooutliers$POROSITY, xlim=c(-4,4), main="GoM_proc_nooutliers", col="red",
     xlab="POROSITY")

# Ensuring the same distribution of Target Variable ORF in Test and Train Data sets
par(mfrow=c(2,1))
hist(GoMTrain_processed$ORF, main="Training ORF distribution",
     xlab="ORF", col="green")

hist(GoMTest_processed$ORF, main="Test ORF distribution",
     xlab="ORF", col="green")

par(mfrow=c(1,1))


####################################################################################
                            #Non Singular Data set
####################################################################################
GoM_processed_Nonsingular<- GoM_processed
GoM_processed_Nonsingular[,c("FSTRU.K","PLAY_TYPE.P1",
                             "CHRONOZONE2.PLIO_UPR","DRIVE.WTR","RES_TYPE.U")]<-NULL


GoMTrain_processed_Nonsingular<-GoM_processed_Nonsingular[train,]

GoMTest_processed_Nonsingular<-GoM_processed_Nonsingular[-train,]

#####################################################################################
#                         MODEL 1: Linear Regression                                #
#####################################################################################

#Removing hightly correlated predictors
 

# Fitting a simple linear regression with no cross validation

Simple_linear_model<-lm(ORF~., data=GoMTrain_processed)
summary(Simple_linear_model)
# Prediction Using Simple Linear Regression
s<-predict(Simple_linear_model, newdata=GoMTest_processed)

# RMSE
RMSE (s, GoMTest_processed$ORF)
# 0.09

# MAE
MAE (s, GoMTest_processed$ORF)
# 0.073

plot( GoMTest_processed$ORF,s, xlab="Original", ylab="Predicted", col="red",
      main="ORF prdicted vs original", xlim=c(0,1), ylim=c(-0.3,1))
abline(0,1, col="blue", lwd=2)


smry<-summary(Simple_linear_model)   # Residual standard error: 0.076, Adjusted R-squared:  0.81
AIC(Simple_linear_model)     # -3674.6
plot(Simple_linear_model)

Simple_linear_model_worst_predictions<-data.frame(Original= GoMTest_processed$ORF,
                                            Predicted= predict(Simple_linear_model, GoMTest_processed))

Simple_linear_model_worst_predictions$Error<-abs(Simple_linear_model_worst_predictions$Original-
                                                   Simple_linear_model_worst_predictions$Predicted)

Worst_predictions<-Simple_linear_model_worst_predictions[order(-Simple_linear_model_worst_predictions$Error),]
head(Worst_predictions)

Best_Predictions<-Simple_linear_model_worst_predictions[order(Simple_linear_model_worst_predictions$Error),]
head(Best_Predictions)

hist(Simple_linear_model_worst_predictions$Error, main="Distribution of absolute error", xlab="Error",
        col="red",breaks=20)

TotalPredictions<-data.frame(Original= GoM_processed$ORF,
                    Simple_linear_model= predict(Simple_linear_model, GoM_processed))



#####################################################################################
#                         MODEL 1.2: Linear Regression  with non linear terms
#####################################################################################

#Removing hightly correlated predictors


# Fitting a simple linear regression with no cross validation



linear_model_inter<-lm(ORF~(.)^2, data=GoMTrain_processed)
sss<-summary(linear_model_inter)
rownames(sss$coefficients)[sss$coefficients[,4]<0.05]
 
summary(linear_model_inter)
# Prediction Using Simple Linear Regression
s<-predict(linear_model_inter, newdata=GoMTest_processed)

# RMSE
RMSE (s, GoMTest_processed$ORF)
#  All interactions 0.09

# MAE
MAE (s, GoMTest_processed$ORF)
#  All interactions0.073

plot( GoMTest_processed$ORF,s, xlab="Original", ylab="Predicted", col="red",
      main="ORF prdicted vs originalm with interactions between Predictors", xlim=c(0,1), ylim=c(-0.3,1))
abline(0,1, col="blue", lwd=2)

TotalPredictions$Simple_with_nonlinear <-  predict(linear_model_inter, GoM_processed)


# MLR only with stat significant predictors-----------------------------------------------

linear_model_limited_inter<-lm(ORF~.+FSTRU.A:Np+ FSTRU.B:Np + 
                        FSTRU.C:Np+ FSTRU.D:Np  + FSTRU.E:Np+
                        FSTRU.K:Np+ FSTRU.A:Mobility_Ratio_endpoint + 
                          FSTRU.B:Mobility_Ratio_endpoint + 
                         FSTRU.C:Mobility_Ratio_endpoint+ 
                          FSTRU.D:Mobility_Ratio_endpoint +  
                         FSTRU.E:Mobility_Ratio_endpoint+
                      FSTRU.K:Mobility_Ratio_endpoint+ 
                        FSTRU.E:CHRONOZONE2.MIOCENE + 
                        FSTRU.E:CHRONOZONE2.PLEISTOCENE + 
                        FSTRU.D:PLAY_TYPE.A1,
                     data=GoMTrain_processed)


limited_inter<-summary(linear_model_limited_inter)
rownames(limited_inter$coefficients)[sss$coefficients[,4]<0.05]

summary(limited_inter)
# Prediction Using Simple Linear Regression
s<-predict(linear_model_limited_inter, newdata=GoMTest_processed)

# RMSE
RMSE (s, GoMTest_processed$ORF)
#  All interactions 0.09

# MAE
MAE (s, GoMTest_processed$ORF)
#  All interactions0.073

plot( GoMTest_processed$ORF,s, xlab="Original", ylab="Predicted", col="red",
      main="ORF prdicted vs originalm with interactions between Predictors", xlim=c(0,1), ylim=c(-0.3,1))
abline(0,1, col="blue", lwd=2)

TotalPredictions$Simple_with_nonlinear <-  predict(linear_model_inter, GoM_processed)
#####################################################################################
#                         MODEL 2: Robust Regression                                #
#####################################################################################


# Fitting a Robust linear regression with no cross validation

Robust_linear_model<-rlm(ORF~., GoMTrain_processed_Nonsingular, psi=psi.bisquare)
summary(Robust_linear_model)
# Prediction Using Simple Linear Regression
s_Robust<-predict(Robust_linear_model, newdata=GoMTest_processed_Nonsingular)

# RMSE
RMSE (s_Robust, GoMTest_processed_Nonsingular$ORF)
# 0.093

# MAE
MAE (s_Robust, GoMTest_processed_Nonsingular$ORF)
# 0.073

plot( GoMTest_processed_Nonsingular$ORF,s_Robust, xlab="Original", ylab="Predicted", col="red",
      main="ORF prdicted vs original", xlim=c(0,1), ylim=c(-0.3,1))
abline(0,1, col="blue", lwd=2)



Robust_linear_model_worst_predictions<-data.frame(Original= GoMTest_processed_Nonsingular$ORF,
                                      Predicted= predict(Robust_linear_model, 
                                                         GoMTest_processed_Nonsingular))

Robust_linear_model_worst_predictions$Error<-abs(Robust_linear_model_worst_predictions$Original-
                                                   Robust_linear_model_worst_predictions$Predicted)

Worst_predictions<-Robust_linear_model_worst_predictions[order(-Robust_linear_model_worst_predictions$Error),]
head(Worst_predictions)

Best_Predictions<-Robust_linear_model_worst_predictions[order(Robust_linear_model_worst_predictions$Error),]
head(Best_Predictions)

hist(Robust_linear_model_worst_predictions$Error, main="Distribution of absolute error", xlab="Error",
     col="red",breaks=20)

TotalPredictions$Robust_linear_model<- predict(Robust_linear_model, GoM_processed)

  
#####################################################################################
#                       Model-3: LASSO                                                      #
#####################################################################################

set.seed(1)
# Running Lasso with 10 folds cross validation
#------------------------------------------------------------------------------------
names(GoMTrain_processed)
lasso_model <- cv.glmnet(as.matrix(GoMTrain_processed[,-31]), as.numeric(GoMTrain_processed$ORF), alpha=1,
                         nfolds = 3,type.measure = "mae")

plot(lasso_model) # Plot of MSE vs log Lambda
plot(lasso_model$glmnet.fit, xvar="lambda", label=TRUE, main="Coefficient Paths")


#1SE model
coef(lasso_model, s=lasso_model$lambda.1se)
Se_num<-which(lasso_model$lambda== lasso_model$lambda.1se)
lasso_model$cvm[Se_num]

# Model with 11 predictors-----------------------------------------------------
num_11<- which(lasso_model$nzero==11)
coef(lasso_model, s=lasso_model$lambda[num_11])
lasso_model$cvm[num_11] 
lasso_model$lambda[num_11]


s_lasso_11<-predict(lasso_model, newx = as.matrix(GoMTest_processed[,-31]), 
                    s=lasso_model$lambda[num_11])

RMSE_lasso_11 <- RMSE(s_lasso_11[,1], GoMTest_processed$ORF)
RMSE_lasso_11

plot( GoMTest_processed$ORF,s_lasso_11[,1], xlab="Original", ylab="Predicted", col="green",
      main="ORF prdicted vs original (LASSO - 11 predictors)", xlim=c(0,1), ylim=c(-0.3,1))

abline(0,1, col="blue", lwd=2)

# Model with 7 predictors-----------------------------------------------------
num_7<- which(lasso_model$nzero==7)
coef(lasso_model, s=lasso_model$lambda[num_7])
lasso_model$cvm[num_7] 
lasso_model$lambda[num_7]

s_lasso_7<-predict(lasso_model, newx = as.matrix(GoMTest_processed[,-31]), 
                    s=lasso_model$lambda[num_7])

RMSE_lasso_7 <- RMSE(s_lasso_7, GoMTest_processed$ORF)
RMSE_lasso_7

plot( GoMTest_processed$ORF,s_lasso_7, xlab="Original", ylab="Predicted", col="green",
      main="ORF prdicted vs original (LASSO - 7 predictors)", xlim=c(0,1), ylim=c(-0.3,1))

abline(0,1, col="blue", lwd=2)


##################################################################################
#        MODEL 4 K Nearest neighbours
##################################################################################
kNNprediction<- data.frame()

for(i in 1:nrow(GoMTest_processed)){
  nN<-nn2(GoMTrain_processed, GoMTest_processed[i,], k=25)
  temp<-mean(GoMTrain_processed[nN$nn.idx,]$ORF)
  kNNprediction[i,1]<-i
  kNNprediction[i,2]<-temp
}

names(kNNprediction)<-c("Index", "ORF")

# RMSE
RMSE (kNNprediction$ORF, GoMTest_processed$ORF)

# MAE
MAE (kNNprediction$ORF, GoMTest_processed$ORF)


plot( GoMTest_processed$ORF,kNNprediction$ORF, xlab="Original", ylab="Predicted", col="purple",
      main="ORF prdicted vs original (kNN-25)", xlim=c(0,1), ylim=c(-0.3,1))
abline(0,1, col="blue", lwd=2)

# k    RMSE    MAE
# 1   0.168    0.129
# 5   0.134    0.108
#10   0.128    0.105
#25   0.127    0.106
#50   0.129    0.108
#100  0.132    0.111
#200  0.135    0.113

kNNprediction$original<- GoMTest_processed$ORF
kNNprediction$Error<-abs(kNNprediction$ORF-kNNprediction$original)
rownames(kNNprediction)<-rownames(GoMTest_processed)

hist(kNNprediction$Error, main="Distribution of absolute error (kNN-25)", xlab="Error",
     col="purple",breaks=20)

#Best predictions
head(kNNprediction[order(kNNprediction$Error),])

#Worst Predictions
head(kNNprediction[order(-kNNprediction$Error),])

#Prediction for all reservoir instances using kNN

for(i in 1:nrow(GoM_processed)){
  nN<-nn2(GoM_processed[-i,], GoM_processed[i,], k=25)
  temp<-mean(GoM_processed[nN$nn.idx,]$ORF)
  TotalPredictions$kNN[i]<-temp
}


#####################################################################################
#                              MODEL 5: Decision Tree                              #
#####################################################################################

# Decision tree based on original data 
DT_model<-rpart(ORF~.,data=GoMTrain_original,                   
                parms=list(split="information"),   
                control=rpart.control(cp=0.003),xval=2)  

plotcp(DT_model)
# cp=0.02 and tree size = 9
par(mfrow=c(1,1))

#Importance fo predictors
barplot(DT_model$variable.importance, main="Variable Importance", col="brown" )
varImp(DT_model)

#Pruning back the tree based on 1 SE
DT_model<-rpart(ORF~.,data=GoMTrain_original,              
                parms=list(split="information"),   
                control=rpart.control(cp=0.01),xval=3) 
# Plot DT
fancyRpartPlot(DT_model)

# Using party and partykit packages
GoMparty<-as.party(DT_model)
plot(GoMparty)



s_DT<-predict(DT_model, newdata=GoMTest_original)

RMSE<-RMSE(s_DT, GoMTest_original$ORF)
RMSE
# 0.13

MAE<-MAE(s_DT, GoMTest_original$ORF)
MAE
# 0.10

# Regression tree of size 24 with Cp 0.005
DT_model_25<-rpart(ORF~.,data=GoMTrain_original,              
                parms=list(split="information"),   
                control=rpart.control(cp=0.005),xval=3) 

s_DT_25<-predict(DT_model_25, newdata=GoMTest_original)


plot( GoMTest_original$ORF,s_DT, xlab="Original", ylab="Predicted", col="brown",
      main="ORF predicted vs original (Regression Tree size-9)", xlim=c(0,0.8), ylim=c(0,0.8))
abline(0,1, col="blue", lwd=2)

#Evaluation of predictions using Regression tree

DT_predictions<-data.frame(Original= GoMTest_original$ORF,
                                            Predicted= predict(DT_model, GoMTest_original))

DT_predictions$Error<-abs(DT_predictions$Original-DT_predictions$Predicted)

Worst_predictions<-DT_predictions[order(-DT_predictions$Error),]
head(Worst_predictions)

Best_Predictions<-DT_predictions[order(DT_predictions$Error),]
head(Best_Predictions)

hist(DT_predictions$Error, main="Distribution of absolute error (Decision Tree)", xlab="Error",
        col="brown",breaks=20)

TotalPredictions$DT_model<-predict(DT_model, GoM_original)

#####################################################################################
#                              MODEL 8: Random Forest                              #
#####################################################################################

RF_model<-randomForest(ORF~., data=GoMTrain_original, ntree=1000, mtry=3)

p_RF<-predict(RF_model, newdata=GoMTest_original)

RMSE<-RMSE(p_RF, GoMTest_original$ORF)
RMSE
# 0.11


MAE<-MAE(p_RF, GoMTest_original$ORF)
MAE
# 0.09
plot( GoMTest_original$ORF,p_RF, xlab="Original", ylab="Predicted", col="blue",
      main="ORF predicted vs original (Random Forest)", xlim=c(0,0.8), ylim=c(0,0.8))
abline(0,1, col="red", lwd=2)

#Evaluation of predictions using Random forest

RF_predictions<-data.frame(Original= GoMTest_original$ORF,
                           Predicted= predict(RF_model, GoMTest_original))

RF_predictions$Error<-abs(RF_predictions$Original-RF_predictions$Predicted)

Worst_predictions<-RF_predictions[order(-RF_predictions$Error),]
head(Worst_predictions)

Best_Predictions<-RF_predictions[order(RF_predictions$Error),]
head(Best_Predictions)

hist(RF_predictions$Error, main="Distribution of absolute error (Random Forest)", xlab="Error",
     col="blue",breaks=20)

TotalPredictions$RF_model<-predict(RF_model, GoM_original)
#####################################################################################
#                              MODEL 6: Neural Network                              #
#####################################################################################

NN_model<-nnet(ORF~., data=GoMTrain_processed, size=2, maxit = 1000)

plot(NN_model)
p<-predict(NN_model, newdata=GoMTest_processed)

RMSE<-RMSE(p, GoMTest_processed$ORF)
RMSE
# 0.084

MAE<-MAE(p, GoMTest_processed$ORF)
MAE
# 0.06
plot( GoMTest_processed$ORF,p, xlab="Original", ylab="Predicted", col="green", pch=16,
      main="ORF predicted vs original (ANN)", xlim=c(0,0.8), ylim=c(0,0.8))
abline(0,1, col="blue", lwd=2)



#Evaluation of predictions using ANN

NN_predictions<-data.frame(Original= GoMTest_original$ORF,
                           Predicted= predict(NN_model, GoMTest_processed))

NN_predictions$Error<-abs(NN_predictions$Original-NN_predictions$Predicted)

Worst_predictions<-NN_predictions[order(-NN_predictions$Error),]
head(Worst_predictions)

Best_Predictions<-NN_predictions[order(NN_predictions$Error),]
head(Best_Predictions)

hist(NN_predictions$Error, main="Distribution of absolute error (ANN)", xlab="Error",
     col="green",breaks=20)

plot(NN_predictions$Predicted, NN_predictions$Original-NN_predictions$Predicted,
     xlab="Fitter values", ylab= "Residuals", main="Residuals Vs Fitted Values")

abline(h=0, lwd=2, col="blue")

TotalPredictions$NN_model<-predict(NN_model, GoM_processed)




################################################################################
#             Mean of all the models
################################################################################
TotalPredictions$mean_prediction= (TotalPredictions$Simple_linear_model+
                                     TotalPredictions$RF_model+TotalPredictions$NN_model)/3

TotalTest<-TotalPredictions[-train,]
RMSE<-RMSE(TotalTest$mean_prediction, GoMTest_processed$ORF)
RMSE
# 0.086

MAE<-MAE(TotalTest$mean_prediction, GoMTest_processed$ORF)
MAE
#0.0687




################################################################################
#                           Ensemble modelling
################################################################################

joint<- TotalPredictions[,c(2,6,7)]
Ensemble_ful<-cbind(GoM_processed, joint)



# Setting seed for reproducibility
set.seed(22)

# Data is split by keeping same distribution of ORF in train and test data
train<-createDataPartition(GoM_original$ORF, p=0.80,list=FALSE)
#Ensemble train and test datasets
GoMTrain_Ensemble<-Ensemble_ful[train,]

GoMTest_Ensemble<-Ensemble_ful[-train,]

# Decision tree based on original data 
#Ensemble_model<- randomForest(ORF~., data=GoMTrain_Ensemble, ntree=1000, mtry=3)

Ensemble_model<-rpart(ORF~.,data=GoMTrain_Ensemble,              
                parms=list(split="information"),   
                control=rpart.control(cp=0.005),xval=3) 


plotcp(Ensemble_model)
# cp=0.02 and tree size = 9
par(mfrow=c(1,1))


# Plot DT
fancyRpartPlot(Ensemble_model)

# Using party and partykit packages
GoMparty<-as.party(Ensemble_model)
plot(GoMparty)


 s_Ensemble<-predict(Ensemble_model, newdata=GoMTest_Ensemble)

 RMSE<-RMSE(s_Ensemble, GoMTest_Ensemble$ORF)
 RMSE
 # 0.0627
 
 MAE<-MAE(s_Ensemble, GoMTest_Ensemble$ORF)
 MAE
 # 0.0465

Ensemble_predictions<-data.frame(Original= GoMTest_Ensemble$ORF,
                            Predicted= predict(Ensemble_model, GoMTest_Ensemble))
 
Ensemble_predictions$Error<-abs(Ensemble_predictions$Original-Ensemble_predictions$Predicted)
 
 Worst_predictions<-Ensemble_predictions[order(-Ensemble_predictions$Error),]
 head(Worst_predictions)
 
 Best_Predictions<-Ensemble_predictions[order(Ensemble_predictions$Error),]
 head(Best_Predictions)
 
 hist(Ensemble_predictions$Error, main="Distribution of absolute error (Ensemble)", xlab="Error",
      col="orange",breaks=20)
 
TotalPredictions$Ensemble<-predict(Ensemble_model, Ensemble_ful)
TotalPredictions$Ensemble_error<-TotalPredictions$Original-TotalPredictions$Ensemble

############################################################################################################

rownames(oilResClean)<-oilResClean$SAND_NAME





################################################################################
#                          Profile generator RF model  Bhcomp
################################################################################
profile_bhcomp<-function(reservoir=1, BHCOMPmin=1, BHCOMPmax=101){
  
  profile<-data.frame(No.of_wells= double(), ORF=double())
  BHCOMPIncrement<-(BHCOMPmax-BHCOMPmin)/10
  ii=1
  
  for(BHCOMP in seq(from=BHCOMPmin, to=BHCOMPmax, by=BHCOMPIncrement)){
    
    temp_res<-GoMTest_original[reservoir,]
    temp_res$BHCOMP<-BHCOMP
    
    name=rownames(GoMTest_processed[reservoir,])
    k=oilResClean[name,]$PERMEABILITY
    area=oilResClean[name,]$OAREA
    dev= k * BHCOMP/area
    temp_res$Dev_factor<-dev
    
    p<-predict(RF_model, temp_res)
    profile[ii,1]<-BHCOMP
    profile[ii,2]<-p
    #print(BHCOMP)
    
    ii=ii+1
  }
  
  print(profile)
  plot(profile$No.of_wells, profile$ORF, ylab="RF", xlab="No of wells", 
       main=paste("Sand name: ", name), pch=16, col="red", type="b" )
  return(profile)
  
}

test<-profile_bhcomp(reservoir= 150, BHCOMPmin=10, BHCOMPmax=80)
#50
#5
#100
#150




################################################################################
#                          Profile generator RF model  Porosity
################################################################################
profile_porosity<-function(reservoir=1, Pormin=0.05, Pormax=0.35){
  
  profile<-data.frame(Porosity= double(), ORF=double())
  PorIncrement<-(Pormax-Pormin)/10
  ii=1
  
  for(Porosity in seq(from=Pormin, to=Pormax, by=PorIncrement)){
    
    temp_res<-GoMTest_original[reservoir,]
    temp_res$POROSITY<-Porosity
    name=rownames(temp_res)
  
    p<-predict(RF_model, temp_res)
    profile[ii,1]<-Porosity
    profile[ii,2]<-p
    #print(BHCOMP)
    
    ii=ii+1
  }
  
  print(profile)
  plot(profile$Porosity, profile$ORF, ylab="RF", xlab="Porosity", 
       main=paste("Sand name: ", name), pch=16, col="blue", type="b" )
  return(profile)
  
}

test<-profile_porosity(reservoir= 40)
test<-profile_porosity(reservoir= 47)
test<-profile_porosity(reservoir= 150)
test<-profile_porosity(reservoir= 106)


################################################################################
#                          Profile generator RF model  permeabiltiy
################################################################################
profile_perm<-function(reservoir=1, Kmin=1, Kmax=1000){
  
  profile<-data.frame(K= double(), ORF=double())
  KIncrement<-(Kmax-Kmin)/10
  ii=1
  
  for(K in seq(from=Kmin, to=Kmax, by=KIncrement)){
    
    temp_res<-GoMTest_original[reservoir,]

    
    name=rownames(GoMTest_processed[reservoir,])
    k=oilResClean[name,]$PERMEABILITY
    
    temp_res$Dev_factor<-temp_res$Dev_factor*K/k
    
    #temp_res$Ng<-temp_res$Ng*K/k
    
    #temp_res$Npc<-temp_res$Npc*sqrt(K)/sqrt(k)
    
    p<-predict(RF_model, temp_res)
    profile[ii,1]<-K
    profile[ii,2]<-p
    #print(BHCOMP)
    
    ii=ii+1
  }
  
  print(profile)
  plot(profile$K, profile$ORF, ylab="RF", xlab="Permeability (md)", 
       main=paste("Sand name: ", name), pch=16, col="brown", type="b" )
  return(profile)
  
}

test<-profile_perm(reservoir= 120, Kmin=10, Kmax=1000)
test<-profile_perm(reservoir= 1, Kmin=10, Kmax=1000)
test<-profile_perm(reservoir= 50, Kmin=10, Kmax=1000)
test<-profile_perm(reservoir= 70, Kmin=10, Kmax=1000)
#50
#5
#100
#150


#Distribution of predictors
par(mfrow=c(3,3))

hist(rCleanSelect$POROSITY, main=NULL, xlab="Porosity")
hist(rCleanSelect$SW, main=NULL, xlab="Sw")
hist(rCleanSelect$BHCOMP, main=NULL, xlab="Wells")
hist(rCleanSelect$Mobility_Ratio_endpoint, main=NULL, xlab="Mob. Ratio")
hist(rCleanSelect$Nalpha, main=NULL, xlab="Nalpha")
hist(rCleanSelect$Npc, main=NULL, xlab="Npc")
hist(rCleanSelect$Dev_factor,main=NULL, xlab="Dev factor")
hist(rCleanSelect$Hetro_factor, main=NULL, xlab="Hetro factor")
hist(rCleanSelect$ORF, main=NULL, xlab="ORF")

par(mfrow=c(1,1))

# Skewness of each predictor
skewness(rCleanSelect$POROSITY)
skewness(rCleanSelect$SW)
skewness(oilResClean$PERMEABILITY)
skewness(rCleanSelect$BHCOMP)
skewness(rCleanSelect$Npc)
skewness(rCleanSelect$Ng)
skewness(rCleanSelect$Nalpha)
skewness(rCleanSelect$Dev_factor)
skewness(rCleanSelect$Hetro_factor)



#Afte transformation
par(mfrow=c(3,2))

hist(GoM_original$BHCOMP, main="Original", xlab="No of Wells", col="blue")
hist(GoM_processed$BHCOMP, main="Transformaed", xlab="No of Wells", col="blue")


hist(GoM_original$Dev_factor, main=NULL, xlab="Dev Factor", col="red")
hist(GoM_processed$Dev_factor, main=NULL, xlab="Dev Factor", col="red")



hist(GoM_original$Npc, main=NULL, xlab="Npc", col="green")
hist(GoM_processed$Npc, main=NULL, xlab="Npc", col="green")


#Distribution of Error
TotalPredictions$Simple_linear_model_Error<-TotalPredictions$Original-TotalPredictions$Simple_linear_model
TotalPredictions$Robust_linear_model_Errror<-TotalPredictions$Original-TotalPredictions$Robust_linear_model
TotalPredictions$kNN_Error<-TotalPredictions$Original-TotalPredictions$kNN
TotalPredictions$DT_model_Error<-TotalPredictions$Original-TotalPredictions$DT_model
TotalPredictions$RF_model_Error<-TotalPredictions$Original-TotalPredictions$RF_model
TotalPredictions$NN_model_Error<-TotalPredictions$Original-TotalPredictions$NN_model
TotalPredictions$Liner_with_inter_error <- TotalPredictions$Original-TotalPredictions$Simple_with_nonlinear

plot(TotalPredictions$Original, TotalPredictions$Simple_linear_model_Error, col=rgb(1,0,0,0.2),
     pch=16, ylab="Error", xlab="ORF-Original", xlim=c(0,1), main="Error vs Original ORF")
points(TotalPredictions$Original, TotalPredictions$kNN_Error, col=rgb(1,1,0,0.2), pch=16)
points(TotalPredictions$Original, TotalPredictions$RF_model_Error, col=rgb(0,0,1,0.2),pch=16)
points(TotalPredictions$Original, TotalPredictions$NN_model_Error, col=rgb(0,1,0,0.2),pch=16)

abline(h=0, lwd=2, col="blue")

legend(0.7,-0.05, c("Multi linear regression","kNN", "Random Forest", "ANN"),
       pch=c(16,16,16,16),col=c("red","yellow","blue","green"))


#Error
d<-density(abs(TotalPredictions$Simple_linear_model_Error))
plot(d, lwd=2, col="red", ylim=c(0,15), main="Distribution of Error", xlab="Error")
dd<-density(abs(TotalPredictions$kNN_Error))
points(dd, lwd=2, col="yellow", type="l")
ddd<-density(abs(TotalPredictions$RF_model_Error))
points(ddd, lwd=2, col="blue", type="l")
dddd<-density(abs(TotalPredictions$NN_model_Error))
points(dddd, lwd=2, col="green", type="l")
ddddd<-density(abs(TotalPredictions$Ensemble_error))
points(dddd, lwd=2, col="orange", type="l")

legend(0.45,15, c("Multi linear regression","kNN", "Random Forest", "ANN"),
       lwd=c(2,2,2,2),col=c("red","yellow","blue","green"))


# Error vs Original plots----------------------------------------------------------------
# Multiple liner regression
plot(GoMTest_processed$ORF, predict(Simple_linear_model, newdata= GoMTest_processed), 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=17  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF - predict(Simple_linear_model, newdata= GoMTest_processed)), 
     col="blue", type="p", 
     pch=2)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(17, 2), col=c("red", "blue"), cex=1)

#Mutiple linear regression with all possible in Interactions------------------------------
mod=linear_model_inter
predPoint<- 15
errorPoint<- 0

plot(GoMTest_processed$ORF, predict(mod, newdata= GoMTest_processed), 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF - predict(mod, newdata= GoMTest_processed)), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


#Random Forest------------------------------
mod=RF_model
predPoint<- 16
errorPoint<- 1

plot(GoMTest_processed$ORF, predict(mod, newdata= GoMTest_original), 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF - predict(mod, newdata= GoMTest_original)), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


#  ANN------------------------------
mod=NN_model
predPoint<- 18
errorPoint<- 

plot(GoMTest_processed$ORF, predict(mod, newdata= GoMTest_processed), 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF - predict(mod, newdata= GoMTest_processed)), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


#  kNN-----------------------------

predPoint<- 19
errorPoint<- 6
  
plot(GoMTest_processed$ORF, kNNprediction$ORF, 
       ylab="Predicted ORF/Error", xlab="Original ORF", 
       #main="Prediction Quality (MLR)",
       col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, kNNprediction$Error, 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


#  LASSO 11 Predictors------------------------------

predPoint<- 16
errorPoint<- 3
  
  plot(GoMTest_processed$ORF, s_lasso_11[,1], 
       ylab="Predicted ORF/Error", xlab="Original ORF", 
       #main="Prediction Quality (MLR)",
       col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF - s_lasso_11[,1]), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


#  Regression tree of size 9 -----------------------------

predPoint<- 16
errorPoint<- 4

plot(GoMTest_processed$ORF, s_DT, 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF-s_DT), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


#  Regression tree of size 24 -----------------------------

predPoint<- 16
errorPoint<- 4

plot(GoMTest_processed$ORF, s_DT_25, 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF-s_DT_25), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))

#  Ensemble model--------------------------------------------------------------


predPoint<- 15
errorPoint<- 11

plot(GoMTest_Ensemble$ORF, s_Ensemble, 
     ylab="Predicted ORF/Error", xlab="Original ORF", 
     #main="Prediction Quality (MLR)",
     col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_Ensemble$ORF, (GoMTest_Ensemble$ORF-s_Ensemble), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))


hist(GoM_original$Nalpha)


#  MLR with non linear interactions-------------------------------------------------------------


mod=linear_model_limited_inter
predPoint<- 15
errorPoint<- 0
  
  plot(GoMTest_processed$ORF, predict(mod, newdata= GoMTest_processed), 
       ylab="Predicted ORF/Error", xlab="Original ORF", 
       #main="Prediction Quality (MLR)",
       col="red", type="p", xlim=c(0,0.7), ylim=c(-0.2,0.8), pch=predPoint  )
abline(a=0, b=1, col="black", lwd=2)
abline(h=0, col="black", lwd=2)

points(GoMTest_processed$ORF, (GoMTest_processed$ORF - predict(mod, newdata= GoMTest_processed)), 
       col="blue", type="p", 
       pch=errorPoint)

legend(0,0.8, c("Prediction", "Error"), 
       pch=c(predPoint, errorPoint), col=c("red", "blue"))