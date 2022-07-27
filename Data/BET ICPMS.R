#Set working directory, load biostats.r, load libraries
setwd('C:/Users/zwhitener/Box/Kerr Lab/Fisheries Science Lab/PERC project/BET and Albacore/ICPMS data')

#save plots
graphics.off()
windows(width=50,  height=26, record=TRUE)
.SavedPlots <- NULL

### Required packages ####
library(plyr)
library(ggplot2)
library(histogram)
library(lattice)
library(MASS)
library(raster)
library(gridExtra)
library(Rmisc)

####Read in data####
rawdata<-(read.csv("ICPMS master_060622.csv")) 
summary(rawdata)

#Remove 'bad data'
dataR<- subset(rawdata, Cut.Out=='0') # removed data deemed as "bad" during ICPMS run
dataS<- subset(dataR, Plus.=='0')

summary(dataS)
data<- dataS

datafish<-subset(data, Type=='Sample')
summary(datafish)

#how to remove "dead" data from samples only?
data<- subset(dataS, X7Li > 0 & X25Mg > 0 & X43Ca > 0 & X55Mn > 0 & X63Cu > 0 & 
              X68Zn > 0 & X88Sr > 0 &  X138Ba > 0  ) #removed values below detection limit
summary(data)



####Data Exploration and Visualization#### 
ggplot(data=datafish, aes(x=Seg, y=X88Sr, group=ID, colour=Type)) +
  geom_line() +
  geom_point()

#Checking individual ID
Checkdata<-subset(data, ID=='2308-2018')
#qplot(Seg,X88Sr, geom='smooth', span =0.5, data=Checkdata,group=Age,colour="blue")+ylim(0,5)
A<-ggplot(data=Checkdata, aes(x=Age.of.Sample, y=X88Sr, group=ID, colour=Age)) +
  geom_line() +
  geom_point()
B<-ggplot(data=Checkdata, aes(x=Age.of.Sample, y=X138Ba, group=ID, colour=Age)) +
  geom_line() +
  geom_point()
C<-ggplot(data=Checkdata, aes(x=Age.of.Sample, y=X55Mn, group=ID, colour=Age)) +
  geom_line() +
  geom_point()
D<-ggplot(data=Checkdata, aes(x=Age.of.Sample, y=X25Mg, group=ID, colour=Age)) +
  geom_line() +
  geom_point()
grid.arrange(A,B,C,D) 


##
#Summary by location
Adata<-subset(datafish, Region=='HAT')
Bdata<-subset(datafish, Region=='NHAT')
Cdata<-subset(datafish, Region=='NMAB')
Ddata<-subset(datafish, Region=='SMAB')



#Point data plot (unsmoothed)

ggplot(data=Adata, aes(x=Seg, y=X88Sr, group=ID, colour=Year.Class)) +
  geom_line() +
  geom_point()

ggplot(data=Ddata, aes(x=Dist.from.core, y=Ba138.Ca48, group=sample_id, colour=Location.Spawning)) +
  geom_line() +
  geom_point()

####STOP 063022##

#Smoothed plot Dist from Core
qplot(Dist.from.core,Sr88.Ca48, geom='smooth', span =0.5, data=Adata,group=sample_id,colour="blue")+ylim(0,5)
#Smoothed plot Age of growth
qplot(Age.of.Growth,Sr88.Ca48, geom='smooth', data=Cdata,group=sample_id,colour="blue")+ylim(0,5)
#newdata<- data[-c(sample_id=='061512_1')]

#Smoothed element data over time based on spawning time
Wdata<-subset(data, Spawning.time=='W')
qplot(Dist.from.core,Sr88.Ca48, geom='smooth', span =0.5, data=Wdata,group=sample_id,colour="blue")+ylim(0,5)
Sdata<-subset(data, Spawning.time=='S')
qplot(Dist.from.core,Sr88.Ca48, geom='smooth', span =0.5, data=Sdata,group=sample_id,colour=Spawning.time)+ylim(0,5)
qplot(Dist.from.core,Sr88.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Spawning.time)+ylim(0,5)

qplot(Dist.from.core,Sr88.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)
qplot(Dist.from.core,Ba138.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)
qplot(Dist.from.core,Mn55.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)
qplot(Dist.from.core,Mg25.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)
qplot(Dist.from.core,Cu63.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)
qplot(Dist.from.core,Zn68.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)
qplot(Dist.from.core,Cd114.Ca48, geom='smooth', span =0.5, data=data,group=sample_id,colour=Location.Spawning)

#####CORE OTOLITH ANALYSIS####
datacore<- subset(data, Age.of.Growth<2 & Age.of.Growth>-2)
summary(datacore)
data_coremed<- ddply(datacore, c("sample_id", "Sex","Spawning.time", "Age", "YEAR","Length", "Weight",
                                 "Capture.Location","Location.Spawning"), summarise,
                     medianMg25=median(Mg25.Ca48,na.rm=TRUE),
                     medianMn55=median(Mn55.Ca48,na.rm=TRUE),
                     medianCu63=median(Cu63.Ca48,na.rm=TRUE),  
                     medianBa138=median(Ba138.Ca48,na.rm=TRUE),
                     medianSr88=median(Sr88.Ca48,na.rm=TRUE),
                     medianZn68=median(Zn68.Ca48,na.rm=TRUE),  
                     medianSr86=median(Sr86.Ca48,na.rm=TRUE),
                     medianSr87=median(Sr87.Ca48,na.rm=TRUE),
                     medianBa137=median(Ba137.Ca48,na.rm=TRUE),
                     CVMg25=cv(Mg25.Ca48,na.rm=TRUE),
                     CVMn55=cv(Mn55.Ca48,na.rm=TRUE),
                     CVCu63=cv(Cu63.Ca48,na.rm=TRUE),  
                     CVZn68=cv(Zn68.Ca48,na.rm=TRUE),  
                     CVBa138=cv(Ba138.Ca48,na.rm=TRUE),
                     CVSr88=cv(Sr88.Ca48,na.rm=TRUE),
                     CVSr86=cv(Sr86.Ca48,na.rm=TRUE),
                     CVSr87=cv(Sr87.Ca48,na.rm=TRUE),
                     CVBa137=cv(Ba137.Ca48,na.rm=TRUE),
                     meanMg25=mean(Mg25.Ca48,na.rm=TRUE),
                     meanMn55=mean(Mn55.Ca48,na.rm=TRUE),
                     meanCu63=mean(Cu63.Ca48,na.rm=TRUE),  
                     meanZn68=mean(Zn68.Ca48,na.rm=TRUE),  
                     meanBa138=mean(Ba138.Ca48,na.rm=TRUE),
                     meanSr88=mean(Sr88.Ca48,na.rm=TRUE),
                     meanSr87=mean(Sr87.Ca48,na.rm=TRUE),
                     meanSr86=mean(Sr86.Ca48,na.rm=TRUE),
                     meanBa137=mean(Ba137.Ca48,na.rm=TRUE))
summary (data_coremed)
data_coremed <- data_coremed[-c(59,73), ]
row.names (data_coremed)=NULL

#Length at age
#Dage<-subset(data_coremed, Age=='5')
#Dage<-subset(data_coremed, Age<'6')
#Len<- aov(Length ~ Age*Spawning.time,data = Dage)
##summary(Len)
#ggplot(data=Dage, aes(x=Age, y=Length, group=Spawning.time, colour=Spawning.time)) +
#geom_point()+ geom_smooth(aes(x=Age, y=Length, group=Spawning.time, colour=Spawning.time))

# Histogram Overlapping Distribution of 
c<-ggplot(subset(data_coremed,!is.na(medianBa138))) + geom_histogram(data=data_coremed,aes(x=medianBa138,fill=Spawning.time),alpha=.5)
d<-ggplot(subset(data_coremed,!is.na(medianSr88))) + geom_histogram(data=data_coremed,aes(x=medianSr88,fill=Spawning.time),alpha=.5)
e<-ggplot(subset(data_coremed,!is.na(medianMg25))) + geom_histogram(data=data_coremed,aes(x=medianMg25,fill=Spawning.time),alpha=.5)
f<-ggplot(subset(data_coremed,!is.na(medianMn55))) + geom_histogram(data=data_coremed,aes(x=medianMg25,fill=Spawning.time),alpha=.5)
grid.arrange(c,d,e,f, ncol=2)

##Oto core element boxplots
library(gridExtra)
j<-ggplot(SRdata_coremed,mapping=aes(y=medianSr88,x=Spawning.time,colour=Capture.Location)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Sr:Ca",x="Spawning Time")+ theme(legend.position="none")
k<-ggplot(data_coremed,mapping=aes(y=medianBa138,x=Spawning.time, colour=Capture.Location)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Ba:Ca",x="Spawning Time")+ theme(legend.position="none")
l<-ggplot(MGdata_coremed,mapping=aes(y=medianMg25,x=Spawning.time, colour=Capture.Location)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Mg:Ca",x="Spawning Time")+ theme(legend.position =c(0.8,0.8))
m<-ggplot(MNdata_coremed,mapping=aes(y=medianMn55,x=Spawning.time, colour=Capture.Location)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Mn:Ca",x="Spawning Time")+ theme(legend.position="none")
a<-ggplot(CUdata_coremed,mapping=aes(y=medianCu63,x=Spawning.time, colour=Capture.Location)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Cu:Ca",x="Spawning Time")+ theme(legend.position="none")
b<-ggplot(ZNdata_coremed,mapping=aes(y=medianZn68,x=Spawning.time, colour=Capture.Location)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Zn:Ca",x="Spawning Time")+ theme(legend.position="none")
grid.arrange(j,k,l,m,a,b, ncol=3)


#Oto core element By Year
ggplot(subset(data_coremed,!is.na(medianSr88))) + geom_smooth(aes(y=medianSr88,x=YEAR, colour=Spawning.time))+geom_point(aes(y=medianSr88,x=YEAR,colour=Spawning.time))
ggplot(subset(data_coremed,!is.na(medianBa138))) + geom_smooth(aes(y=medianBa138,x=YEAR, colour=Spawning.time))+geom_point(aes(y=medianBa138,x=YEAR,colour=Spawning.time))
ggplot(subset(data_coremed,!is.na(medianMg25))) + geom_smooth(aes(y=medianMg25,x=YEAR, colour=Spawning.time))+geom_point(aes(y=medianMg25,x=YEAR,colour=Spawning.time))
ggplot(subset(data_coremed,!is.na(medianMn55))) + geom_smooth(aes(y=medianMn55,x=YEAR, colour=Spawning.time))+geom_point(aes(y=medianMn55,x=YEAR,colour=Spawning.time))

# Kernal density contours for the baseline data 
ggplot(data_coremed) + geom_density2d(aes(y=medianMg25, x=medianBa138, col=Spawning.time)) + geom_point(aes(y=medianMg25, x=medianBa138, col=Spawning.time))
ggplot(data_coremed) + geom_density2d(aes(y=medianSr88, x=medianBa138, col=Spawning.time)) + geom_point(aes(y=medianSr88, x=medianBa138, col=Spawning.time))
ggplot(data_coremed) + geom_density2d(aes(y=medianSr88, x=medianMn55, col=Spawning.time)) + geom_point(aes(y=medianSr88, x=medianMn55, col=Spawning.time))

#histograms
par(mfcol = c(2, 6))
for (k in 9:14) {
  j0 <- names(data_coremed)[k]
  br0 <- seq(min(data_coremed[, k]), max(data_coremed[, k]), le = 11)
  x0 <- seq(min(data_coremed[, k]), max(data_coremed[, k]), le = 50)
  for (i in 1:2) {
    i0 <- levels(data_coremed$Spawning.time)[i]
    x <- data_coremed[data_coremed$Spawning.time == i0, j0]
    hist(x, br = br0, proba = T, col = grey(0.8), main = i0,
         xlab = j0)
    lines(x0, dnorm(x0, mean(x), sd(x)), col = "red", lwd = 2)
  }
}
#Bivariate Scatterplots
library(ade4)
par(mar = c(0, 0, 0, 0))
pan1 <- function(x, y, ...) {
  xy <- cbind.data.frame(x, y)
  s.class(xy, data_coremed$Spawning.time, include.ori = F, add.p = T, clab = 1.5,
          col = c("red", "blue"), cpoi = 2, csta = 0.5)
}
pairs(data_coremed[, 9:14], panel = pan1)

####Core oto MULTIVARIATE ANALYSIS####
####MANOVA
library(MASS)
library(MVN)
library(mvnormtest)
library(mvoutlier)
library(nlme)
#MV <- manova(cbind(log(medianMg25),medianMn55,log(medianCu63),log(medianZn68),medianSr88,log(medianBa138))~Capture.Location*Spawning.time,  data=data_coremed)
#MV <- manova(cbind(CVMg25,CVMn55,CVCu63,CVZn68,CVSr88,CVBa138)~Capture.Location*Spawning.time,  data=data_coremed)
MV <- manova(cbind(log(medianMg25),log(medianMn55),log(medianCu63),medianSr88,log(medianBa138))~Capture.Location*Spawning.time*YEAR,  data=data_coremed)
summary(MV, test="Pillai")
plot(MV)
summary.aov(MV)
MV$outliers
aq.plot(data_coremed[, 9:14])
x<-data_coremed[, 9:14]
y<-na.omit(x)
res <- mvoutlier.CoDa(y)
plot(res,onlyout=TRUE,bw=FALSE,which="parallel",symb=TRUE,symbtxt=TRUE,transp=0.3)

# Are data normally distributed?-Univariate Normality
ggplot(data=data.frame(qqnorm( data_coremed$meanSr88, plot=F),Stock=data_coremed$Spawning.time), mapping=aes(x=x,y=y)) +geom_point()+geom_smooth(method="lm",se=FALSE) + facet_wrap(~Stock) 
ggplot(data=data.frame(qqnorm( data_coremed$medianBa138, plot=F),Stock=data_coremed$Spawning.time), mapping=aes(x=x,y=y)) +geom_point()+geom_smooth(method="lm",se=FALSE) + facet_wrap(~Stock) 
ggplot(data=data.frame(qqnorm( data_coremed$medianMg25, plot=F),Stock=data_coremed$Spawning.time), mapping=aes(x=x,y=y)) +geom_point()+geom_smooth(method="lm",se=FALSE) + facet_wrap(~Stock) 
ggplot(data=data.frame(qqnorm( data_coremed$medianMn55, plot=F),Stock=data_coremed$Spawning.time), mapping=aes(x=x,y=y)) +geom_point()+geom_smooth(method="lm",se=FALSE) + facet_wrap(~Stock) 
ggplot(data=data.frame(qqnorm( data_coremed$medianCu63, plot=F),Stock=data_coremed$Spawning.time), mapping=aes(x=x,y=y)) +geom_point()+geom_smooth(method="lm",se=FALSE) + facet_wrap(~Stock) 

##Multivariate normality
MVdata<-data_coremed
MVdata[, 8:12]<-log(data_coremed[, 8:12]+1)
#MVdata2<-na.omit(logMVdata) #ommit NAs
mardiaTest <- mardiaTest(MVdata[, 8:13], qqplot = FALSE) #multivariate normality 
mardiaTest
hzTest <- hzTest(MVdata[, 8:12], qqplot = TRUE) #multivariate normality
hzTest
# Adjusted Mahalanobis distance
outlier <- mvOutlier(MVdata, qqplot = TRUE, method = "adj.quan")
outlier
#MVdata3 <- MVdata2[-c(49,76,50,39,42), ]
uniPlot(MVdata, type = "qqplot") # creates univariate Q-Q plots
uniPlot(MVdata, type = "histogram") # creates univariate histograms
uniNorm(MVdata, type = "SW", desc = TRUE) #performs univariate normality test

#library(biotools)
#boxM(data_coremed, Spawning.time)
###random Forest Classification

### Random Forest Classification Model####
data_coremed2=na.omit(data_coremed)
summary(data_coremed2)
set.seed(50)
base.rf <- randomForest(Spawning.time~medianMg25+medianMn55+medianCu63+medianZn68+medianSr88+medianBa138+CVBa138+CVSr88+CVMg25+CVMn55+CVCu63+CVZn68, data=data_coremed,importance=TRUE,
                        proximity=TRUE,oob.prox=T)
print(base.rf)
imp <- importance(base.rf)

#revised Model
set.seed(50)
rf2 <- randomForest(Spawning.time~medianMg25+medianMn55+medianSr88+medianBa138+CVBa138+CVSr88+CVMg25+CVCu63, data=data_coremed,importance=TRUE,
                    proximity=TRUE,oob.prox=T)
print(rf2)
imp <- importance(rf2)

summary(data_coremed)
set.seed(50)
rf2 <- randomForest(Spawning.time~medianMg25+medianMn55+medianSr88+medianBa138+medianCu63, data=data_coremed,importance=TRUE,
                    proximity=TRUE,oob.prox=T)
print(rf2)
imp <- importance(rf2)


####Core oto Discriminant Function Analysis####
##PCA
library(ade4)
par(mfcol = c(3, 4))
pca1 <- dudi.pca(data_coremed[, 10:15], scannf = FALSE)

#Linear discriminant analysis
dis1 <- discrimin(pca1, data_coremed$Spawning.time, scannf = FALSE)
names(dis1)
dis1
plot(dis1)

#Histogram plots of discrimination based on all parameters
library(MASS)
dis2 <- lda(as.matrix(log(data_coremed[, 10:14])), data_coremed$Spawning.time, color=Spawning.time)
names(dis2)
dis2
plot(dis2)

dis3 <- lda(as.matrix(log(data_coremed[, 10:14])), data_coremed$Capture.Location)
names(dis3)
dis3
plot(dis3)

dis4 <- lda(as.matrix(log(data_coremed[, 10:14])), data_coremed$Location.Spawning)
names(dis4)
dis4
plot(dis4)

#Manual exploration of discrimination by  each var
y<-(data_coremed[ ,10:15])
x <- log(data_coremed[ ,10:13])
z <- cbind(x,data_coremed$medianSr88)
#y<-subset(data_coremed,select=c("medianSr88", "medianBa138","medianMg25"))
grp<-(data_coremed[ ,3])
grp2<-(data_coremed[ ,8])
grp3<-(data_coremed[ ,9])
y.mat<-as.matrix(y)
summary(lm(y.mat~grp))

#Stepwise DA
out<-greedy.wilks(z,grp)
out2<-greedy.wilks(z,grp2)
out3<-greedy.wilks(z,grp3)

#DFA1 <- lda(Spawning.time~meanMg25+meanMn55+meanCu63+meanZn68+meanSr88+meanBa138, data=data_coremed,CV=TRUE)
#DFA2 <- lda(Spawning.time~medianMg25+medianMn55+medianCu63+medianZn68+medianSr88+medianBa138, data=data_coremed,CV=TRUE)
#plot(DFA1)

#DFA Spawning Group with Jacknife leave-one-out
#DFA <- lda(Spawning.time~meanSr88+meanBa138, data=data_coremed,CV=TRUE)
#DFA <- lda(Spawning.time~meanMg25+meanMn55+meanCu63+meanZn68+meanSr88+meanBa138, data=data_coremed,CV=TRUE)
#DFA <- lda(Spawning.time~medianMg25+medianMn55+log(medianCu63)+medianZn68+medianSr88+medianBa138, data=data_coremed,CV=TRUE)
#DFA <- lda(Spawning.time~log(medianMg25)+log(medianMn55)+log(medianZn68)+log(medianSr88)+log(medianBa138)+CVSr88+CVBa138, data=data_coremed,CV=TRUE)
#DFA <- lda(Spawning.time~medianMg25+medianMn55+medianZn68+medianSr88+log(medianBa138), data=data_coremed,CV=TRUE)
#DFA <- lda(Spawning.time~(medianMg25)+log(medianMn55)+(medianSr88)+log(medianBa138), data=data_coremed, CV=TRUE)

#Final Spawning Time DFA model
DFA <- lda(Spawning.time~medianSr88+log(medianBa138)+log(medianMg25)+log(medianMn55),data=data_coremed,CV=TRUE)
summary(DFA)
#lda.values<- predict(lda, data_coremed[ ,10:15])$class
ct <- table(data_coremed$Spawning.time,DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))
plot(DFA)


#Final Capture Location DFA model
DFA <- lda(Capture.Location~medianSr88+log(medianMn55)+log(medianCu63),data=data_coremed,CV=TRUE)
summary(DFA)
lda.values<- predict(lda, data_coremed[ ,10:15])$class
ct <- table(data_coremed$Spawning.time,DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

#scores<- prcomp(data_coremed[ ,9:14])
#DFA Capture Location
#DFA <- lda(Capture.Location~medianMg25+medianMn55+medianZn68+medianSr88+medianBa138, data=data_coremed,CV=TRUE)
DFA <- lda(Capture.Location~medianMg25+medianMn55+medianZn68+medianSr88+medianBa138+CVSr88+CVBa138, data=data_coremed,CV=TRUE)
summary(DFA)
lda.values<- predict(lda, data_coremed[ ,10:15])$class
ct <- table(data_coremed$Capture.Location, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

#Bivariate plot
plot(data_coremed$medianMg25, data_coremed$medianSr88, col=Spawning.time)
text(data_coremed$medianMg25, data_coremed$medianSr88,data_coremed$Spawning.time, cex=0.7, pos=4, col="red")


#plot(DFA, dimen=1, type="both") # fit from lda

#plot(DFA, panel = panel.lda, cex = 0.7, dimen=1,
abbrev = FALSE, xlab = "LD1", ylab = "LD2")

# DFA with location*spawning
DFA <- lda(Location.Spawning~log(medianMg25)+log(medianMn55)+medianSr88+log(medianBa138)+log(medianCu63), data=data_coremed,CV=TRUE)
ct <- table(data_coremed$Location.Spawning, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))


pairs(data_coremed[ ,10:15])
DFA$scaling

DFA.pred <- predict(DFA)
DFA.pred
plot(DFA)

### CORE UNIVARIATE ANALYSIS ####
data_coremed
####CORE COPPER
CUdata_coremed <- data_coremed[-c(27,60,121), ]
coreCu63<- aov(log(medianCu63)  ~ Capture.Location*Spawning.time,data =CUdata_coremed)
summary(coreCu63)
plot(coreCu63)
#Test for Equal Variance
bartlett.test(log(medianCu63)~Spawning.time, data=CUdata_coremed)
#boxplot
bwplot(medianCu63 ~ Spawning.time|Capture.Location , data=CUdata_coremed )
ggplot(CUdata_coremed,mapping=aes(y=medianCu63,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()
# summary
summarySE(CUdata_coremed, measurevar="medianCu63", groupvars=c("Spawning.time", "Location.Spawning"))

####CORE ZINC
ZNdata_coremed <- data_coremed[-c(2,3,39, 42,100, 112), ]
coreZn68<- aov((medianZn68)  ~ Capture.Location*Spawning.time,data = ZNdata_coremed)
summary(coreZn68)
plot(coreZn68)
#Test for Equal Variance
bartlett.test(medianZn68~Spawning.time, data=ZNdata_coremed)
# test for Normality
shapiro.test(ZNdata_coremed$medianZn68)
#boxplot
bwplot(medianZn68 ~ Spawning.time|Capture.Location , data=data_coremed )
# summary
summarySE(ZNdata_coremed, measurevar="medianZn68", groupvars=c("Spawning.time", "Location.Spawning"))

####CORE MAGNESIUM
MGdata_coremed<- data_coremed[-c(47,96,103), ]
coreMg25<- aov(log(medianMg25)  ~ Capture.Location*Spawning.time,data =MGdata_coremed)
summary(coreMg25)
plot(coreMg25)
#Test for Equal Variance
bartlett.test(log(medianMg25)~Spawning.time, data=MGdata_coremed)
#boxplot
ggplot(MGdata_coremed,mapping=aes(y=medianMg25,x=Spawning.time, colour=Spawning.time)) + geom_boxplot()
# summary
summarySE(MGdata_coremed, measurevar="medianMg25", groupvars=c("Spawning.time", "Location.Spawning"))

####CORE MANGENESE
MNdata_coremed <- data_coremed[-c(20,31,53), ]
coreMn55<- aov(log(medianMn55)  ~ Capture.Location*Spawning.time,data = MNdata_coremed)
summary(coreMn55)
plot(coreMn55)
#Test for Equal Variance
bartlett.test(medianMn55~Spawning.time, data=MNdata_coremed)
#boxplot
ggplot(MNdata_coremed,mapping=aes(y=medianMn55,x=Spawning.time, colour=Spawning.time)) + geom_boxplot()
# summary
summarySE(MNdata_coremed, measurevar="medianMn55", groupvars=c("Spawning.time", "Location.Spawning"))

####CORE STRONTIUM
SRdata_coremed<- data_coremed[-c(16,78), ]# removed outlier
coreSr88<- aov(medianSr88  ~ Capture.Location*Spawning.time,data =SRdata_coremed)
summary(coreSr88)
plot(coreSr88)
#Test for Equal Variance
bartlett.test(medianSr88~Spawning.time, data=SRdata_coremed)
#boxplot
ggplot(SRdata_coremed,mapping=aes(y=medianSr88,x=Location.Spawning, colour=Spawning.time)) + geom_boxplot()
# summary
summarySE(SRdata_coremed, measurevar="medianSr88", groupvars=c("Spawning.time", "Location.Spawning"))

####CORE BARIUM 
#logtransform 
coreBa138<- aov(log(medianBa138)  ~ Capture.Location*Spawning.time,data = data_coremed)
summary(coreBa138)
plot(coreBa138)
#Test for Equal Variance
bartlett.test(log(medianBa138)~Spawning.time, data=data_coremed)
#Box plot
ggplot(data_coremed,mapping=aes(y=medianBa138,x=Spawning.time, colour=Spawning.time)) + geom_boxplot()
# summary
summarySE(data_coremed, measurevar="medianBa138", groupvars=c("Spawning.time", "Location.Spawning"))

#CV in Sr:Ca
SRcvdata_coremed <- data_coremed[-c(120), ]
CVSr88<- aov((CVSr88)  ~ Capture.Location*Spawning.time,data =SRcvdata_coremed)
summary(CVSr88)
plot(CVSr88)
#Test for Equal Variance
bartlett.test((CVSr88)~Spawning.time, data=SRcvdata_coremed)

#CV in Ba:Ca
BAcvdata_coremed <- data_coremed[-c(36), ]
CVBa138<- aov((CVBa138)  ~ Capture.Location*Spawning.time,data =BAcvdata_coremed)
summary(CVBa138)
plot(CVBa138)
#Test for Equal Variance
bartlett.test((CVBa138)~Spawning.time, data=BAcvdata_coremed)

#CV in Mg:Ca
MGcvdata_coremed <- data_coremed[-c(112), ]
CVMg25<- aov(log(CVMg25)  ~ Capture.Location*Spawning.time,data =MGcvdata_coremed)
summary(CVMg25)
plot(CVMg25)
#Test for Equal Variance
bartlett.test((CVMg25)~Spawning.time, data=MGcvdata_coremed)

#CV in Mn:Ca
MNcvdata_coremed <- data_coremed[-c(10), ]
CVMn55<- aov(log(CVMn55)  ~ Capture.Location*Spawning.time,data =MNcvdata_coremed)
summary(CVMn55)
plot(CVMn55)
#Test for Equal Variance
bartlett.test((CVMn55)~Spawning.time, data=MNcvdata_coremed)

#CV in Cu:Ca
CUcvdata_coremed <- data_coremed[-c(20,129), ]
CVCu63<- aov(log(CVCu63)  ~ Capture.Location*Spawning.time,data =CUcvdata_coremed)
summary(CVCu63)
plot(CVCu63)
#Test for Equal Variance
bartlett.test((CVMn55)~Spawning.time, data=MNcvdata_coremed)

##summary graphs
library(gridExtra)
w<-bwplot(medianSr88 ~ Spawning.time, data=SRdata_coremed)
x<-bwplot (medianMg25 ~ Spawning.time , data=MGdata_coremed)
y<-bwplot(medianMn55 ~ Spawning.time , data=MNdata_coremed)
z<-bwplot(medianBa138 ~ Spawning.time , data=BAdata_coremed)
grid.arrange(w, x, y,z, ncol=2)

##summary graphs
library(gridExtra)
w<-bwplot(medianSr88 ~ Location.Spawning , data=data_coremed)
x<-bwplot(log(medianMg25) ~ Location.Spawning , data=data_coremed)
y<-bwplot(log(medianMn55) ~ Location.Spawning , data=data_coremed)
z<-bwplot(medianBa138 ~ Location.Spawning , data=data_coremed)
grid.arrange(w, x, y,z, ncol=2)
#####PLOTS FOR CORE UNIVARIATE RESPONSE#####
library (ggplot2)
#medianSr88
Sr<-ggplot(data_coremed,aes(Spawning.time,medianSr88,fill=factor(Spawning.time)))+ 
  geom_boxplot(size=1)+ylim(c(1.5, 3.5))+ylab("Sr:Ca (mmole/mole)") +xlab("Spawning Group")+
  theme_classic()+ theme(axis.line = element_line(colour = "black"),
                         axis.text=element_text(size=16),axis.title=element_text(size=18, face="bold"),legend.position="none")
#medianBa138
Ba<-ggplot(data_coremed,aes(Spawning.time,medianBa138,fill=factor(Spawning.time)))+ 
  geom_boxplot(size=1)+ylab("Ba:Ca (mmole/mole)") +xlab("Spawning Group")+
  theme_classic()+ theme(axis.line = element_line(colour = "black"),
                         axis.text=element_text(size=16),axis.title=element_text(size=18, face="bold"),legend.position="none")
#medianMg25
Mg<-ggplot(data_coremed,aes(Spawning.time,medianMg25,fill=factor(Spawning.time)))+ 
  geom_boxplot(size=1)+ylab("Mg:Ca(mmole/mole)") +xlab("Spawning Group")+
  theme_classic()+ theme(axis.line = element_line(colour = "black"),
                         axis.text=element_text(size=16),axis.title=element_text(size=18, face="bold"),legend.position="none")+
  scale_fill_manual(breaks = c("S", "W"),values=c( "blue", "green"))

#medianMn
Mn<-ggplot(data_coremed,aes(Spawning.time,medianMn55,fill=factor(Spawning.time)))+ 
  geom_boxplot(size=1)+ylab("Mn:Ca (mmole/mole)") +xlab("Spawning Group")+
  theme_classic()+ theme(axis.line = element_line(colour = "black"),
                         axis.text=element_text(size=16),axis.title=element_text(size=18, face="bold"),legend.position="none")
grid.arrange(Sr, Ba, Mg,Mn, ncol=2)

##### Temperature signal in Core Oto data#####
#Overall trend
library (Rmisc)
temp_all <- summarySE(data_coremed, measurevar="medianSr88", groupvars=c("YEAR"))

qplot(YEAR,medianSr88, geom='smooth', span =0.5, data=data_coremed, xlab="Year", ylab="Sr:Ca", xlim=c(2008,2011))+ theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.75),
        axis.line.y = element_line(color="black", size = 0.75),cex.lab=1.5)

data_coremed2<-subset(data_coremed, )
Sr_year<- aov(log(medianSr88)~YEAR*Spawning.time,data = data_coremed)
summary(Sr_year)
plot(Sr_year)

#trend by origin
qplot(YEAR,medianSr88, geom='smooth', span =0.5, data=data_coremed,group=Spawning.time,colour=Spawning.time, xlab="Year", ylab="Sr:Ca", xlim=c(2008,2011) )+ theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.75),
        axis.line.y = element_line(color="black", size = 0.75))

qplot( YEAR,medianSr88, geom='smooth', span =0.5, data=data_coremed,group=Spawning.time,colour=Spawning.time, xlim=c(2008,2011), xlab="Year", ylab= "SrCa")
qplot( YEAR,medianBa138, geom='smooth', span =0.5, data=data_coremed,group=Spawning.time,colour=Spawning.time, xlim=c(2008,2011), xlab="Year", ylab= "BaCa")
qplot( YEAR,medianMg25, geom='smooth', span =0.5, data=data_coremed,group=Spawning.time,colour=Spawning.time, xlim=c(2008,2011), xlab="Year", ylab= "MgCa")
qplot( YEAR,medianMn55, geom='smooth', span =0.5, data=data_coremed,group=Spawning.time,colour=Spawning.time, xlim=c(2008,2011), xlab="Year", ylab= "MnCa")

ggplot(data_coremed,aes(x=YEAR, y=medianSr88,group=YEAR))+
  ##Boxplot across years grouped by origin
  ggplot(data_coremed,aes(x=YEAR, y=medianSr88,group=YEAR))+
  geom_boxplot(aes(group = interaction(factor(YEAR), Spawning.time),
                   fill = Spawning.time))+  
  theme_classic()+ylab("Sr:Ca(mmol/mole)")+xlab("Year")

ggplot(data_coremed,aes(x=YEAR, y=medianSr88,group=YEAR))+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 2), size = 1)
#Mg Plot across years
ggplot(data_coremed,aes(x=YEAR, y=medianMg25,group=YEAR))+
  geom_boxplot(aes(group = interaction(factor(YEAR), Spawning.time),
                   fill = Spawning.time))
#####WHOLE OTOLITH ANALYSIS####
data2<- subset(data, Age.of.Growth>-2)
summary(data2)

data_whole<- ddply(data2, c("sample_id", "Sex","Spawning.time", "Age", "Capture.Location", "Location.Spawning"), summarise,
                   medianMg25=(median(Mg25.Ca48)),
                   medianMn55=(median(Mn55.Ca48)),
                   medianBa138=(median (Ba138.Ca48)),
                   medianSr88=(median(Sr88.Ca48)),
                   medianZn68=(median(Zn68.Ca48)), 
                   CVBa138=cv(Ba137.Ca48),
                   medianSr86=(median(Sr86.Ca48)),
                   medianSr87=(median(Sr87.Ca48)),
                   medianCd114=(median(Cd114.Ca48)),  
                   medianBa137=(median(Ba137.Ca48)),
                   medianCu63=(median(Cu63.Ca48)), 
                   CVMg25=cv(Mg25.Ca48),
                   CVMn55=cv(Mn55.Ca48),
                   CVCu63=cv(Cu63.Ca48),  
                   CVZn68=cv(Zn68.Ca48),  
                   CVSr88=cv(Sr88.Ca48),
                   CVSr86=cv(Sr86.Ca48),
                   CVSr87=cv(Sr87.Ca48),
                   CVBa137=cv(Ba137.Ca48),
                   meanMg25=mean(Mg25.Ca48),
                   meanMn55=mean(Mn55.Ca48),
                   meanCu63=mean(Cu63.Ca48),  
                   meanZn68=mean(Zn68.Ca48),  
                   meanSr86=mean(Sr86.Ca48),
                   meanSr87=mean(Sr87.Ca48),
                   meanSr88=mean(Sr88.Ca48),
                   meanCd114=mean(Cd114.Ca48),  
                   meanBa137=mean(Ba137.Ca48),
                   meanBa138=mean(Ba138.Ca48))
summary (data_whole)

####WHOLE OTO SR:CA UNIVARIATE#
SRdata_whole <- data_whole[-c(59,72,74), ]
wholeSr88<- aov((medianSr88)  ~ Capture.Location*Spawning.time,data =SRdata_whole)
summary(wholeSr88)
plot(wholeSr88)
#Test for Equal Variance
bartlett.test(medianSr88~Spawning.time, data=SRdata_coremed)

wholeCVSr88<- aov(log(CVSr88)  ~ Capture.Location*Spawning.time,data =SRdata_whole)
summary(wholeCVSr88)
plot(wholeCVSr88)
#Test for Equal Variance
bartlett.test(CVSr88~Spawning.time, data=SRdata_coremed)

####WHOLE OTO BA:CA UNIVARIATE#
BAdata_whole <- data_whole[-c(59,72,126), ]
wholeBa138<- aov(log(medianBa138)  ~ Capture.Location*Spawning.time,data =BAdata_whole)
summary(wholeBa138)
plot(wholeBa138)
#Test for Equal Variance
bartlett.test(medianBa138~Spawning.time, data=BAdata_whole)

BACVdata_whole <- data_whole[-c(58,59,72,90,126), ]
wholeCVBa138<- aov(log(CVBa138) ~ Capture.Location*Spawning.time,data =BACVdata_whole)
summary(wholeCVBa138)
plot(wholeCVBa138)
#Test for Equal Variance
bartlett.test(CVBa138~Spawning.time, data=BACVdata_whole)

####WHOLE OTO MG:CA UNIVARIATE#
MGdata_whole <- data_whole[-c(59,72,97,104), ]
wholeMg25<- aov(log(medianMg25)  ~ Capture.Location*Spawning.time,data =MGdata_whole)
summary(wholeMg25)
plot(wholeCVMg25)
#Test for Equal Variance
bartlett.test(medianMg25~Spawning.time, data=MGdata_whole)

MGCVdata_whole <- data_whole[-c(32,41,58,59, 113), ]
wholeCVMg25<- aov(log(CVMg25)  ~ Capture.Location*Spawning.time,data =MGCVdata_whole)
summary(wholeCVMg25)
plot(wholeCVMg25)
#Test for Equal Variance
bartlett.test(log(CVMg25)~ Spawning.time, data=MGCVdata_whole)

####WHOLE OTO MN:CA UNIVARIATE#
MNdata_whole <- data_whole[-c(59,30,72,98,99), ]
wholeMn55<- aov(log(medianMn55)  ~ Capture.Location*Spawning.time,data =MNdata_whole)
summary(wholeMn55)
plot(wholeMn55)
#Test for Equal Variance
bartlett.test(medianMg25~Spawning.time, data=MNdata_whole)

MNCVdata_whole <- data_whole[-c(10,33,58), ]
wholeCVMn55<- aov(log(CVMn55)  ~ Capture.Location*Spawning.time,data =MNCVdata_whole)
summary(wholeCVMn55)
plot(wholeCVMn55)
#Test for Equal Variance
bartlett.test(log(CVMn55)~Spawning.time, data=MNCVdata_whole)

####WHOLE OTO CU:CA UNIVARIATE#
CUdata_whole <- data_whole[-c(59,72,122), ]
wholeCu63<- aov(log(medianCu63)  ~ Capture.Location*Spawning.time,data =CUdata_whole)
summary(wholeCu63)
plot(wholeCu63)
#Test for Equal Variance
bartlett.test(medianCu63~Spawning.time, data=CUdata_whole)

CUCVdata_whole <- data_whole[-c(58), ]
wholeCVCu63<- aov(log(CVCu63)  ~ Capture.Location*Spawning.time,data =CUCVdata_whole)
summary(wholeCVCu63)
plot(wholeCVCu63)
#Test for Equal Variance
bartlett.test(CVCu63~Spawning.time, data=CUCVdata_whole)

##Oto whole element boxplots
library(gridExtra)
#plotdata_whole <- data_whole[-c(59,72,74), ]
w<-ggplot(SRdata_whole,mapping=aes(y=medianSr88,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Sr:Ca", x="Capture location")+ theme(legend.position="none")
x<-ggplot(BAdata_whole,mapping=aes(y=medianBa138,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Ba:Ca",x="Capture location")+theme(legend.position =c(0.8,0.8))
y<-ggplot(MGdata_whole,mapping=aes(y=medianMg25,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Mg:Ca",x="Capture location")+ theme(legend.position="none")
z<-ggplot(MNdata_whole,mapping=aes(y=medianMn55,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Mn:Ca",x="Capture location")+ theme(legend.position="none")
#a<-ggplot(CUdata_whole,mapping=aes(y=medianCu63,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Cu:Ca",x="Capture location")+ theme(legend.position="none")
#b<-ggplot(Zndata_whole,mapping=aes(y=medianZn68,x=Capture.Location, colour=Spawning.time)) + geom_boxplot()+theme_classic(base_size=14)+labs(y="Zn:Ca",x="Capture location")+ theme(legend.position="none")
grid.arrange(w, x, y,z, ncol=2)

####Whole OTO MANOVA####
MANOVAdata_whole <- data_whole[-c(59,72), ]
#Manova <- manova(cbind(medianMg25,medianMn55,medianCu63,medianSr88,medianBa138)~Capture.Location*Spawning.time,  data=data_whole)
Manova <- manova(cbind(log(medianMg25),log(medianMn55),log(medianCu63),medianSr88,log(medianBa138))~Capture.Location*Spawning.time, data=MANOVAdata_whole)
summary(Manova)
summary.aov(Manova)
#library(biotools)
#boxM(MANOVAdata_whole, Spawning.time)

####CV Whole oto MANOVA####
MANOVAdata_whole <- data_whole[-c(58), ]
Manova <- manova(cbind(log(CVMg25),log(CVMn55),log(CVCu63), log(CVSr88),log(CVBa138))~Capture.Location*Spawning.time,  data=data_whole)
summary(Manova)
summary.aov(Manova)

#histograms
par(mfcol = c(2, 6))
for (k in 8:13) {
  j0 <- names(data_whole)[k]
  br0 <- seq(min(data_whole[, k]), max(data_whole[, k]), le = 11)
  x0 <- seq(min(data_whole[, k]), max(data_whole[, k]), le = 50)
  for (i in 1:2) {
    i0 <- levels(data_whole$Spawning.time)[i]
    x <- data_whole[data_whole$Spawning.time == i0, j0]
    hist(x, br = br0, proba = T, col = grey(0.8), main = i0,
         xlab = j0)
    lines(x0, dnorm(x0, mean(x), sd(x)), col = "red", lwd = 2)
  }
}
####Whole oto Discriminant Analysis####
#Manualexploration of discrimination by  each var
y<-(data_whole[ ,10:15])
x <- log(data_coremed[ ,10:13])
z <- cbind(x,data_coremed$medianSr88)
#y<-subset(data_coremed,select=c("medianSr88", "medianBa138","medianMg25"))
grp<-(data_coremed[ ,3])
grp2<-(data_coremed[ ,8])
grp3<-(data_coremed[ ,9])
y.mat<-as.matrix(y)
summary(lm(y.mat~grp))

##FINAL DFA model
DFA <- lda(Spawning.time~log(medianMg25)+log(medianMn55)+(medianSr88)+log(medianBa138), data=data_whole,CV=TRUE)
ct <- table(data_whole$Spawning.time, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

DFA <- lda(Capture.Location~log(medianCu63), data=data_whole,CV=TRUE)
ct <- table(data_whole$Capture.Location, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

DFA <- lda(Location.Spawning~medianMg25+medianMn55+medianCu63+medianZn68+medianSr88+medianBa138+CVBa138, data=data_whole,CV=TRUE)
DFA <- lda(Location.Spawning~log(medianCu63)+medianSr88, data=data_whole,CV=TRUE)
ct <- table(data_whole$Location.Spawning, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))


#### Random Forest Classification Model: Whole oto
data_whole2=na.omit(data_whole)
set.seed(52)
rf <- randomForest(Spawning.time~medianMg25+medianMn55+medianCu63+medianZn68+medianSr88+medianBa138+CVBa138+CVSr88+CVMg25+CVMn55+CVCu63+CVZn68, data=data_whole2,importance=TRUE,
                   proximity=TRUE,oob.prox=T)
print(rf)
imp <- importance(rf)

#revised Model
set.seed(152)
rf2 <- randomForest(Spawning.time~medianMg25+medianCu63+medianSr88+medianBa138+CVBa138+CVSr88+CVMg25+CVCu63, data=data_whole2,importance=TRUE,
                    proximity=TRUE,oob.prox=T)
print(rf2)
imp <- importance(rf2)

####ANNUAL OTOLITH ANALYSIS####
####DATA INPUT#####
#read in data
data3 <- subset(data, Age.of.Growth>'0.5'& Age.of.Growth<'6' )
datayr<- ddply(data3, c("sample_id","Sex","Spawning.time", "Age", "YEAR",
                        "Capture.Location", "Age.of.Growth", "Location.Spawning"), summarise,
               medianMg25=median(Mg25.Ca48),
               medianMn55=median(Mn55.Ca48),
               medianCu63=median(Cu63.Ca48),  
               medianZn68=median(Zn68.Ca48),  
               medianBa138=median(Ba138.Ca48),
               medianSr88=median(Sr88.Ca48),
               medianSr86=median(Sr86.Ca48),
               medianSr87=median(Sr87.Ca48),
               medianCd114=median(Cd114.Ca48),  
               medianBa137=median(Ba137.Ca48),
               meanMg25=mean(Mg25.Ca48),
               meanMn55=mean(Mn55.Ca48),
               meanCu63=mean(Cu63.Ca48),  
               meanZn68=mean(Zn68.Ca48),  
               meanSr86=mean(Sr86.Ca48),
               meanSr87=mean(Sr87.Ca48),
               meanSr88=mean(Sr88.Ca48),
               meanCd114=mean(Cd114.Ca48),  
               meanBa137=mean(Ba137.Ca48),
               meanBa138=mean (Ba138.Ca48))
summary (datayr)
datayr$Age.of.Growth<- as.factor(datayr$Age.of.Growth)


####AGE Specific ANOVAs by Ratio
#Subset BA  by age for ANOVA
BAdatayr<- subset(datayr, Age.of.Growth=='5')
summary (BAdatayr)
row.names(BAdatayr)=NULL
#BAdatayr2<- BAdatayr[-c(72), ]
#BAdatayr3<- BAdatayr[-c(14,58,72), ]
#BAdatayr4<- BAdatayr[-c(51), ]
BAdatayr5<- BAdatayr[-c(83), ]
#ANOVA
ageBa138<- aov(log(medianBa138) ~ Spawning.time*Capture.Location,data = BAdatayr5)
summary(ageBa138)
plot(ageBa138)

#Subset SR by age for ANOVA
SRdatayr<- subset(datayr, Age.of.Growth=='5')
summary (SRdatayr)
row.names(SRdatayr)=NULL
SRdatayr2<- SRdatayr[-c(72), ]
SRdatayr3<- SRdatayr[-c(58,72,96), ]
SRdatayr4<- SRdatayr[-c(51), ]
SRdatayr5<- SRdatayr[-c(83), ]
#ANOVA
ageSR88<- aov(log(medianSr88) ~ Spawning.time*Capture.Location,data = SRdatayr)
summary(ageSR88)
plot(ageSR88)

#Subset Mn by age for ANOVA
MNdatayr<- subset(datayr, Age.of.Growth=='5')
summary (MNdatayr)
row.names(MNdatayr)=NULL
MNdatayr2<- MNdatayr[-c(28,30, 72), ]
MNdatayr3<- MNdatayr[-c(27,29,26,58,61,72), ]
MNdatayr4<- MNdatayr[-c(51), ]
MNdatayr5<- MNdatayr[-c(83), ]
#ANOVA
ageMN55<- aov(log(medianMn55) ~ Spawning.time*Capture.Location,data = MNdatayr)
summary(ageMN55)
plot(ageMN55)

#Subset MG by age for ANOVA
MGdatayr<- subset(datayr, Age.of.Growth=='5')
summary (MGdatayr)
row.names(MGdatayr)=NULL
MGdatayr2<- MGdatayr[-c(72), ]
MGdatayr3<- MGdatayr[-c(58,72), ]
MGdatayr4<- MGdatayr[-c(51,93), ]
MGdatayr5<- MGdatayr[-c(82), ]
#ANOVA
ageMG25<- aov(log(medianMg25) ~ Spawning.time*Capture.Location,data = MGdatayr4)
summary(ageMG25)
plot(ageMG25)
#################
#Pairwise Comparison
tukey.test <- TukeyHSD(ageBa138)
#Test for Equal Variance
bartlett.test(log(medianBa138)~Location.Spawning, data=BAdatayr)
#Boxplot of Elements across age
ggplot(BAdatayr,aes(x=Age.of.Growth, y=log(medianBa138),group=Age.of.Growth))+
  geom_boxplot(aes(group = interaction(factor(Age.of.Growth), Location.Spawning),
                   fill = Location.Spawning))+theme_classic(base_size=20)+ scale_fill_manual(values=c("red", "blue", "pink", "light blue"))+
  theme(legend.position =c(0.8,0.9))
#Lowess fit 
qplot(Age.of.Growth,medianBa138, geom='smooth', span =0.5, data=BAdatayr,group=Spawning.time,colour=Spawning.time, xlab="Year", ylab= "BaCa")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))
qplot(Age.of.Growth,medianBa138, geom='smooth', span =0.5, data=BAdatayr,group=Location.Spawning,colour=Location.Spawning, xlab="Age", ylab= "BaCa")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))

####AGE Strontium 
SRdatayr<- datayr[-c(274,276,277,337,338,339), ]
#Summarize
Sr <- summarySE(datayr, measurevar="medianSr88", groupvars=c("Age.of.Growth", "Spawning.time"))
#ANOVA
ageSr88<- aov(log(medianSr88)  ~ Age.of.Growth*Spawning.time*Capture.Location,data = SRdatayr)
summary(ageSr88)
#Plots to examine normality, variance, outliers
plot(ageSr88)
#Pariwise comparison
tukey.test <- TukeyHSD(ageSr88)
TukeyHSD(ageSr88, which = 'Spawning.time')
#Boxplot of Elements across age
ggplot(SRdatayr,aes(x=Age.of.Growth, y=log(medianSr88),group=Age.of.Growth))+
  geom_boxplot(aes(group = interaction(factor(Age.of.Growth), Location.Spawning),
                   fill = Location.Spawning))+theme_classic(base_size=20)+ scale_fill_manual(values=c("red", "blue", "pink", "light blue"))+
  theme(legend.position =c(0.8,0.9))
#Lowess fit 
qplot(Age.of.Growth,medianSr88, geom='smooth', span =0.5, data=SRdatayr,group=Spawning.time,colour=Spawning.time, xlab="Age", ylab= "SrCa")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))
qplot(Age.of.Growth,medianSr88, geom='smooth', span =0.5, data=SRdatayr,group=Location.Spawning,colour=Location.Spawning, xlab="Age", ylab= "SrCa")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))

#Test for Equal Variance
bartlett.test(log(medianSr88)~Location.Spawning, data=SRdatayr)

####AGE MAGNESIUM 
MGdatayr<- datayr[-c(277,274,276,337,338,339,502,503), ]
#Summarize
Mg <- summarySE(MGdatayr, measurevar="medianMg25", groupvars=c("Age.of.Growth", "Capture.Location"))
#ANOVA
ageMg25<- aov(log(medianMg25)  ~ Age.of.Growth*Spawning.time*Capture.Location,data = MGdatayr)
summary(ageMg25)
#Plots to examine normality, variance, outliers
plot(ageMg25)
#Pariwise comparison
tukey.test <- TukeyHSD(ageMg25)
#Test for Equal Variance
bartlett.test(log(medianMg25)~Location.Spawning, data=MGdatayr)
#Boxplot of Elements across age
ggplot(MGdatayr,aes(x=Age.of.Growth, y=log(medianMg25),group=Age.of.Growth))+
  geom_boxplot(aes(group = interaction(factor(Age.of.Growth), Location.Spawning),
                   fill = Location.Spawning))+theme_classic(base_size=20)+ scale_fill_manual(values=c("red", "blue", "pink", "light blue"))+
  theme(legend.position =c(0.8,0.9))
#Lowess fit 
qplot(Age.of.Growth,medianMg25, geom='smooth', span =0.5, data=MGdatayr,group=Spawning.time,colour=Spawning.time, xlab="Age", ylab= "Mg25")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))
qplot(Age.of.Growth,medianMg25, geom='smooth', span =0.5, data=MGdatayr,group=Location.Spawning,colour=Location.Spawning, xlab="Age", ylab= "Mg25")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))

####AGE MANGENESE
MNdatayr<- datayr[-c(88,130,131,132, 139,140, 338,339,274,276,277,291,296,337), ]
#Summarize
Mn <- summarySE(MNdatayr, measurevar="medianMn55", groupvars=c("Age.of.Growth", "Spawning.time"))
#ANOVA
ageMn55<- aov(log(medianMn55)  ~ Age.of.Growth*Spawning.time*Capture.Location,data = MNdatayr)
summary(ageMn55)
#Plots to examine normality, variance, outliers
plot(ageMn55)
#Pariwise comparison
TukeyHSD(ageMn55)
#Test for Equal Variance
bartlett.test(log(medianMn55)~Location.Spawning, data=MNdatayr)
#Boxplot of Elements across age
ggplot(MNdatayr,aes(x=Age.of.Growth, y=log(medianMn55),group=Age.of.Growth))+
  geom_boxplot(aes(group = interaction(factor(Age.of.Growth), Location.Spawning),
                   fill = Location.Spawning))+theme_classic(base_size=20)+ scale_fill_manual(values=c("red", "blue", "pink", "light blue"))+
  theme(legend.position =c(0.8,0.9))
#Lowess fit 
qplot(Age.of.Growth,medianMn55, geom='smooth', span =0.5, data=MNdatayr,group=Spawning.time,colour=Spawning.time, xlab="Age", ylab= "Mn55")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))
qplot(Age.of.Growth,medianMn55, geom='smooth', span =0.5, data=MNdatayr,group=Location.Spawning,colour=Location.Spawning, xlab="Age", ylab= "Mn55")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))

####AGE Copper
CUdatayr<- datayr[-c(121,122,210,274,271,272,276,277,337,338,339,398,399,554,555,556,557), ]
#Summarize
Cu <- summarySE(CUdatayr, measurevar="medianCu63", groupvars=c("Age.of.Growth", "Spawning.time"))

ageCu63<- aov(log(medianCu63)  ~ Age.of.Growth*Spawning.time*Capture.Location,data = CUdatayr)
summary(ageCu63)
#Plots to examine normality, variance, outliers
plot(ageCu63)
#Test for Equal Variance
bartlett.test(log(medianCu63)~Location.Spawning, data=CUdatayr)
#Boxplot of Elements across age
ggplot(CUdatayr,aes(x=Age.of.Growth, y=log(medianCu63),group=Age.of.Growth))+
  geom_boxplot(aes(group = interaction(factor(Age.of.Growth), Location.Spawning),
                   fill = Location.Spawning))+theme_classic(base_size=20)+ scale_fill_manual(values=c("red", "blue", "pink", "light blue"))+
  theme(legend.position =c(0.8,0.9))
#Lowess fit 
qplot(Age.of.Growth,medianCu63, geom='smooth', span =0.5, data=CUdatayr,group=Spawning.time,colour=Spawning.time, xlab="Age", ylab= "Cu63")+theme_classic(base_size=20)+theme(legend.position =c(0.8,0.9))



####MANOVA####
library(MASS)
library(mvnormtest)
MV <- manova(cbind(medianMg25,medianMn55,medianCu63,medianZn68,medianSr88,medianBa138)~Capture.Location*Spawning.time*Age.of.Growth,  data=datayr)
summary(MV)
summary.aov(MV)
mshapiro.test(MV)

MV <- manova(cbind(medianMg25,medianMn55,medianCu63,medianZn68,medianSr88,medianBa138)~Age.of.Growth,  data=datayr2)
summary(MV)
#library(biotools)
#boxM(data_coremed, Spawning.time)

####Discriminant Function Analysis####
library(MASS)
DFA1 <- lda(Spawning.time~meanMg25+meanMn55+meanCu63+meanZn68+meanSr88+meanBa138, data=data_coremed,CV=TRUE)
DFA2 <- lda(Capture.Location*Spawning.time*Age.of.Growth~medianMg25+medianMn55+medianCu63+medianZn68+medianSr88+medianBa138, data=datayr,CV=TRUE)

#DFA Spawning Group
DFA <- lda(Spawning.time~medianMg25+medianMn55+medianZn68+medianSr88+medianBa138, data=data_coremed,CV=TRUE)
lda.values <- predict(lda, data_coremed[ ,8:13])$class

ct <- table(data_coremed$Spawning.time, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))


plot(data_coremed$medianMn55, data_coremed$medianSr88)
text(data_coremed$medianMn55, data_coremed$medianSr88,data_coremed$Spawning.time, cex=0.7, pos=4, col="red")


plot(DFA, dimen=1, type="both") # fit from lda

plot(DFA, panel = panel.lda, cex = 0.7, dimen=1,
     abbrev = FALSE, xlab = "LD1", ylab = "LD2")

DFA <- lda(Location.Spawning~medianMg25+medianMn55+medianZn68+medianSr88+medianBa138, data=data_coremed,CV=TRUE)
ct <- table(data_coremed$Location.Spawning, DFA$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))


pairs(data_coremed[ ,8:13])
DFA$scaling

DFA.pred <- predict(DFA)
DFA.pred
plot(DFA)



library(MASS)
lda <- lda(Spawning.time~medianMg25+medianMn55+medianZn68+medianSr88+medianBa138, data=data_coremed)
lda.values <- predict(lda, data_coremed[ ,8:13])$class
lda.values$x[,1] 
class <- predict(lda)$class
ct <- table(lda.values,data_coremed$Spawning.time)

#create a histogram of the discriminant function values
ldahist(data = lda.values$x[,1], g=data_coremed$Spawning.time)
ldahist(data = wine.lda.values$x[,1], g=wine$V1)

#create a scatterplot of the discriminant function values
plot(lda.values$x[,1], type="n", xlim=c(0,30), ylab=c("LDA Axis 1"))
text(lda.values$x[,1], row.names(data_coremed),  col=c(as.numeric(class)+10))
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")

#######NOT SURE
#Outdata<-subset(data, Sr88.Ca48>6) #id outliers in data
#data2<-data[!(data$sample_id =="121214_22"),] #remove outlier
#data3<-data2[!(data$sample_id =="061512_1"),] #remove outlier
# data<- rawdata[!(sample_id='121214_2'& sample_id='060513_7')]
