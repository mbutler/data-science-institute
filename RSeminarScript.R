## ----echo=FALSE,results='hide'-------------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=8, results='markup',cache=TRUE, autodep=TRUE)
knitr::opts_knit$set(progress=FALSE)

## ----eval=FALSE----------------------------------------------------------
## setwd("~/testDirectory/")
## #use getwd() to see working directory
## source("testScript.R")

## ----echo=FALSE----------------------------------------------------------
X=rnorm(10);Y=rnorm(10);DataSet=rnorm(10)

## ------------------------------------------------------------------------
ls()

## ----eval=FALSE----------------------------------------------------------
## save(X,file="test.RData")

## ----eval=FALSE----------------------------------------------------------
## save.image("test.RData")

## ----eval=FALSE----------------------------------------------------------
## load("test.RData")

## ----eval=FALSE----------------------------------------------------------
## getwd()
## setwd("~/newdir/")

## ------------------------------------------------------------------------
#install.packages("lme4")

## ------------------------------------------------------------------------
library(lme4)
#require(lme4)

## ------------------------------------------------------------------------
require(asdf);print("code is still running!")

## ------------------------------------------------------------------------
5L
as.integer(5.0)

## ------------------------------------------------------------------------
5.25
pi
as.numeric(1L)

## ------------------------------------------------------------------------
"hello world"
as.character(1.0)

## ------------------------------------------------------------------------
TRUE
1 > 2
as.logical(1L)

## ------------------------------------------------------------------------
NA
is.na(NA)

## ------------------------------------------------------------------------
NULL
is.null(NULL)

## ----eval=FALSE,echo=FALSE-----------------------------------------------
## rm(list=ls())

## ------------------------------------------------------------------------
x1= 1/3
x2 <- 1/6
1/2 -> x3
class(x1)

## ------------------------------------------------------------------------
(x1+x2-2*x3)*0.1

## ------------------------------------------------------------------------
v1=c(1,2,4,6)
v2=1:4
( v3=0:3*5+1 )

## ------------------------------------------------------------------------
( v1 = c(x1=1,x2=2,x3=4,x4=6) )
names(v2) = c("x1","x2","x3","x4"); print(v2)
names(v3) = paste("x",1:4,sep=""); print(v3)

## ------------------------------------------------------------------------
v1[4]
v1["x4"]
v1[3:4]
v1[c(1,2,4)]
v1[-c(1,2,4)]

## ------------------------------------------------------------------------
( v4 = c("Hello","World") )
class(v4)
( v5 = rep(TRUE,4) )
class(v5)

## ------------------------------------------------------------------------
( fac1= factor(rep(1:3,3),labels=paste("Drug",1:3,sep="")) )

## ------------------------------------------------------------------------
levels(fac1) = paste("Drug",3:1,sep="")
fac1

## ------------------------------------------------------------------------
( mat1 = matrix(0L,2,3) )
( mat2 = matrix(1:6,2,3) )
( mat3 = matrix(1:6,2,3,byrow=TRUE) )

## ------------------------------------------------------------------------
rownames(mat2)=paste("row",1:2)
colnames(mat2)=paste("column",1:3)
mat2

## ------------------------------------------------------------------------
mat2[2,3]
mat2["row 2","column 3"]
mat2[6]
mat2[,c(1,3)]

## ------------------------------------------------------------------------
require("Matrix")
mat4 = matrix(0,1000,50)
mat5 = Matrix(0,1000,50,sparse=TRUE)
class(mat4)
object.size(mat4)
class(mat5)
object.size(mat5)

## ----eval=FALSE----------------------------------------------------------
## require("bigmatrix");require("bigalgebra")
## mat6 = matrix(1.0,5000,5000)
## mat7 = as.big.matrix(mat6)
## object.size(mat6);object.size(mat7)
## #200MB vs. 0.6KB
## system.time(mat6%*%mat6)
## #76.44sec
## system.time(mat7%*%mat7)
## #82.03sec

## ------------------------------------------------------------------------
( arr1 = array(1:12,c(2,3,2)) )

## ------------------------------------------------------------------------
arr1[2,1,1]
arr1[1,,2]

## ----echo=FALSE,results='hide'-------------------------------------------
set.seed(1) 

## ------------------------------------------------------------------------
( df1 = data.frame(id=1:6,
          gender=factor(rep(c("M","F"),3)),
          treatment1=factor(rep(LETTERS[1:3],2)),
          treatment2=factor(rep(LETTERS[1:2],each=3)),
          response=rnorm(6)) )

## ------------------------------------------------------------------------
df1$treatment1
df1[,"treatment1"]
df1[,3]

## ------------------------------------------------------------------------
ls1 = list()
ls1[[1]] = "Hello World"
ls1[[2]] = pi
ls1[[3]] = matrix(1:4,2,2)
print(ls1)

## ------------------------------------------------------------------------
DRfun = function(dose){
  return( 1-exp(-2.18E-04*dose) )
}
DRfun(3.18E+03)
DRfun(1+03)

## ------------------------------------------------------------------------
DRfun = function(dose,K){
  return( 1-exp(-K*dose) )
}
DRfun(1e3,2e-4)

## ------------------------------------------------------------------------
df1[which(df1$response<0),]

## ------------------------------------------------------------------------
with(df1,response[which(treatment1=="C")])

## ------------------------------------------------------------------------
with(df1,by(response,gender,mean))

## ------------------------------------------------------------------------
for(i in 1:4){
  print(i^2)
}

## ------------------------------------------------------------------------
x = c(1,2,4,8,16)
for(i in x){
  print(log2(i))
}

## ------------------------------------------------------------------------
count=0
while(count<5){
  print(count)
  count = count + 1
}

## ------------------------------------------------------------------------
set.seed(1)
x=0
while(x<0.5){
  x = rnorm(1)
  print(x)
}

## ------------------------------------------------------------------------
set.seed(1)
uu = runif(1)
if(uu < 0.5){
  print("heads")
}else{
  print("tails")
}

## ------------------------------------------------------------------------
set.seed(1)
uu = runif(1)
result = ifelse(uu < 0.5, "heads", "tails")
print(result)

## ----eval=FALSE----------------------------------------------------------
## sleepstudy = read.csv("http://myweb.uiowa.edu/dksewell/sleepstudy.csv")
## test = read.table("http://myweb.uiowa.edu/dksewell/sleepstudy.csv",
##                         sep=",",header=TRUE)
## all.equal(sleepstudy,test)

## ------------------------------------------------------------------------
head(sleepstudy)

## ------------------------------------------------------------------------
tail(sleepstudy)

## ----eval=FALSE----------------------------------------------------------
## View(sleepstudy)  #In RStudio only

## ----eval=FALSE----------------------------------------------------------
## print(sleepstudy)

## ------------------------------------------------------------------------
colnames(sleepstudy)

## ------------------------------------------------------------------------
str(sleepstudy)

## ------------------------------------------------------------------------
dim(sleepstudy)

## ------------------------------------------------------------------------
blowdown = alr3::blowdown

## ----eval=FALSE----------------------------------------------------------
## ?blowdown

## ------------------------------------------------------------------------
summary(blowdown)

## ------------------------------------------------------------------------
summary(blowdown$D)

## ------------------------------------------------------------------------
with(blowdown,by(D,y,summary))

## ------------------------------------------------------------------------
colMeans(blowdown[,-4])
#rowMeans()

## ------------------------------------------------------------------------
apply(blowdown[,-4],2,mean)

## ------------------------------------------------------------------------
apply(blowdown[,-4],2,function(x) mean(abs(x-mean(x))))

## ------------------------------------------------------------------------
blowdownScaled = scale(blowdown[,-4])
apply(blowdownScaled,2,
      function(x) round(c(mean=mean(x),sd=sd(x)),10))

## ------------------------------------------------------------------------
require("reshape2")
###Long form:
head(sleepstudy)

## ------------------------------------------------------------------------
###Wide form:
sleepWide = acast(sleepstudy,Subject~factor(Days),
                  value.var="Reaction")
dim(sleepWide)
head(round(sleepWide,1))

## ------------------------------------------------------------------------
###Back to tall form:
sleepTall = melt(sleepWide,value.name="Reaction")
dim(sleepTall)
head(round(sleepTall,1))

## ----eval=F--------------------------------------------------------------
## write.csv(sleepWide,file="sleepWide.csv")
## write.table(sleepWide,file="sleepWide.txt",sep="\t")

## ------------------------------------------------------------------------
X = model.matrix(~Days,data=sleepstudy)
class(X)
dim(X)

## ------------------------------------------------------------------------
qr(X)$rank

## ------------------------------------------------------------------------
XtX = t(X)%*%X
all.equal(XtX,crossprod(X,X))  #Uses less memory

## ------------------------------------------------------------------------
XtXInv = solve(XtX)
# XtXgenInv = MASS::ginv(Xtx)
# XtXInv = chol2inv(chol(XtX))

## ------------------------------------------------------------------------
XtXInv%*%t(X)%*%sleepstudy$Reaction
lm(Reaction~Days,data=sleepstudy)$coef

## ----eval=FALSE----------------------------------------------------------
## ?USJudgeRatings

## ----eval=FALSE----------------------------------------------------------
## judgesCent = scale(as.matrix(USJudgeRatings[,-1]),scale=FALSE)
## Sigma = cov(judgesCent)
## eigs = eigen(Sigma)
## Scores = judgesCent %*% eigs$vectors[,1:2]
## plot(Scores,pch=16,ylim=c(-1,1)*2.5,xlab="",ylab="")
## arrows(0,0,eigs$vec[,1]*3,eigs$vec[,2]*3)
## text(eigs$vec[,1:2]*3,labels=colnames(judgesCent),
##      adj=c(1,0))

## ----height=4,echo=FALSE,results='hide'----------------------------------
judgesCent = scale(as.matrix(USJudgeRatings[,-1]),scale=FALSE)
Sigma = cov(judgesCent)
eigs = eigen(Sigma)
Scores = judgesCent %*% eigs$vectors[,1:2]
plot(Scores,pch=16,ylim=c(-1,1)*2.5,xlab="",ylab="")
arrows(0,0,eigs$vec[,1]*3,eigs$vec[,2]*3)
text(eigs$vec[,1:2]*3,labels=colnames(judgesCent),
     adj=c(1,0))

## ----height=4.5----------------------------------------------------------
yMeans = as.numeric(with(blowdown,by(y,SPP,mean)))
barplot(yMeans,names.arg=levels(blowdown$SPP))

## ----eval=FALSE----------------------------------------------------------
## par(mfrow=c(1,2))
## pie(yMeans,labels=levels(blowdown$SPP),
##     main="Survival rates")
## plotrix::pie3D(yMeans,labels=levels(blowdown$SPP),
##                explode=0.1,
##                col=rgb(20/256,c(1:9*20)/256,120/256))
## par(mfrow=c(1,1))

## ----height=5,echo=FALSE-------------------------------------------------
pie(yMeans,labels=levels(blowdown$SPP))

## ----height=5,echo=FALSE-------------------------------------------------
plotrix::pie3D(yMeans,labels=levels(blowdown$SPP),
               explode=0.1,theta=pi/4,
               col=rgb(20/256,c(1:9*20)/256,120/256))

## ----height=4.5----------------------------------------------------------
boxplot(blowdown$D)

## ----height=4.5----------------------------------------------------------
boxplot(blowdown$D,main="Blowdown data",ylab="Diameter")

## ----height=3,results='hide'---------------------------------------------
par(mfrow=c(1,2))
with(blowdown,by(D,y,boxplot,main="Blowdown data",
                 ylab="Diameter"))
par(mfrow=c(1,1))

## ----height=4------------------------------------------------------------
boxplot(blowdown$D~blowdown$y,main="Blowdown data",
        ylab="Diameter",names=c("survive","died"),
        cex.lab=1.5,cex.axis=1.5)

## ----height=4.5----------------------------------------------------------
hist(blowdown$D,main="Blowdown data",xlab="",
        ylab="Diameter",cex.lab=1.5,cex.axis=1.5)

## ----eval=F--------------------------------------------------------------
## plot(density(blowdown$D),ylab="Diameter",cex.lab=1.5,
##      cex.axis=1.5,main="blowdown data")

## ----height=5,echo=F-----------------------------------------------------
plot(density(blowdown$D),ylab="Diameter",cex.lab=1.5,
     cex.axis=1.5,main="blowdown data")

## ----eval=F--------------------------------------------------------------
## hist(blowdown$D,main="Blowdown data",xlab="",freq=FALSE,
##         ylab="Diameter",cex.lab=1.5,cex.axis=1.5)
## lines(density(blowdown$D,bw=2.5),col="red",lwd=2)

## ----height=5,echo=F-----------------------------------------------------
hist(blowdown$D,main="Blowdown data",xlab="",freq=FALSE,
        ylab="Diameter",cex.lab=1.5,cex.axis=1.5)
lines(density(blowdown$D,bw=2.5),col="red",lwd=2)

## ----height=4.5----------------------------------------------------------
plot(ecdf(blowdown$D),main="Emp. Cum. Distn Func.",
     xlab="",cex.lab=1.5)

## ------------------------------------------------------------------------
DRfun = function(dose){
  return( 1-exp(-2.18E-04*dose) )
}

## ----height=4,echo=FALSE,results='hide'----------------------------------
curve(DRfun(x),from=0,to=50000,xlab="Dose",
      ylab="Probability of Infection",cex.lab=1.5,
      cex.axis=1.5,lwd=2)

## ----eval=FALSE----------------------------------------------------------
## set.seed(1)
## x1 = rnorm(100)
## x2 = rgamma(100,5,10)
## par(mfrow=c(1,2))
## qqnorm(x1);qqline(x1)
## qqnorm(x2);qqline(x2)

## ----height=4.5,results='hide',echo=FALSE--------------------------------
set.seed(123)
x1 = rnorm(100)
x2 = rgamma(100,5,10)
par(mfrow=c(1,2))
qqnorm(x1);qqline(x1)
qqnorm(x2);qqline(x2)

## ----eval=FALSE----------------------------------------------------------
## par(mfrow=c(1,2))
## qqplot(qgamma(ppoints(500),5,10),x1)
## qqline(x1,distribution=function(p)qgamma(p,5,10))
## qqplot(qgamma(ppoints(500),5,10),x2)
## qqline(x2,distribution=function(p)qgamma(p,5,10))

## ----echo=FALSE,results='hide'-------------------------------------------
par(mfrow=c(1,2))
qqplot(qgamma(ppoints(500),5,10),x1)
qqline(x1,distribution=function(p)qgamma(p,5,10))
qqplot(qgamma(ppoints(500),5,10),x2)
qqline(x2,distribution=function(p)qgamma(p,5,10))

## ----eval=FALSE----------------------------------------------------------
## class(ldeaths)
## plot(ldeaths,ylab="# Deaths")
## plot(stl(ldeaths,s.window="periodic"))

## ----echo=FALSE,results='hide'-------------------------------------------
class(ldeaths)
plot(ldeaths,ylab="# Deaths")
# plot(stl(ldeaths,s.window="periodic"))

## ----echo=FALSE,results='hide'-------------------------------------------
class(ldeaths)
# plot(ldeaths,ylab="# Deaths")
plot(stl(ldeaths,s.window="periodic"))

## ----height=4.5----------------------------------------------------------
pairs(airquality)

## ----height=4------------------------------------------------------------
plot(Ozone~Solar.R,data=airquality,xlab="Solar Radiation",
     ylab="Ozone",pch=16)

## ----eval=FALSE----------------------------------------------------------
## CEXs = seq(0.5,5,length.out=500)
## CEXs = with(airquality,
##             CEXs[length(CEXs)*(max(Wind)+1-Wind)/
##                    (max(Wind)+1)])
## with(airquality,
##      plot(Ozone~Solar.R,xlab="Solar Radiation",
##      ylab="Ozone",pch=16,cex=CEXs,col=Month,
##      cex.lab=1.5,cex.axis=1.5))

## ----height=5,echo=FALSE-------------------------------------------------
CEXs = seq(0.5,5,length.out=500)
CEXs = with(airquality,
            CEXs[length(CEXs)*(max(Wind)+1-Wind)/
                   (max(Wind)+1)])
with(airquality,
     plot(Ozone~Solar.R,xlab="Solar Radiation",
     ylab="Ozone",pch=16,cex=CEXs,col=Month,
     cex.lab=1.5,cex.axis=1.5))

## ----eval=FALSE----------------------------------------------------------
## with(airquality,sapply(unique(Month),
##       function(x){ind=which(Month==x);
##                   abline(lm(Ozone[ind]~Solar.R[ind]),
##                          col=x,lwd=3)}))
## with(airquality,legend("topleft",lwd=rep(4,5),cex=1.5,
##        col=unique(Month),legend=unique(Month)))

## ----height=5,echo=FALSE,results='hide'----------------------------------
with(airquality,
     plot(Ozone~Solar.R,xlab="Solar Radiation",
     ylab="Ozone",pch=16,cex=CEXs,col=Month,
     cex.lab=1.5,cex.axis=1.5))
with(airquality,sapply(unique(Month),
      function(x){ind=which(Month==x);
                  abline(lm(Ozone[ind]~Solar.R[ind]),
                         col=x,lwd=3)}))
with(airquality,legend("topleft",lwd=rep(4,5),cex=1.5,
       col=unique(Month),legend=unique(Month)))

## ----eval=FALSE----------------------------------------------------------
## plot(sleepTall$Reaction~sleepTall$Var2,pch=16,
##      xlab="Days of sleep deprivation",
##      ylab="Reaction time",cex.lab=1.5,cex.axis=1.5)
## for(i in 1:nrow(sleepWide)){
##   lines(sleepWide[i,]~c(0:9))
## }
## lines(colMeans(sleepWide)~c(0:9),col="blue",lwd=4,lty=2)

## ----height=4.75,echo=FALSE,results='hide'-------------------------------
plot(sleepTall$Reaction~sleepTall$Var2,pch=16,
     xlab="Days of sleep deprivation",
     ylab="Reaction time",cex.lab=1.5,cex.axis=1.5)
for(i in 1:nrow(sleepWide)){
  lines(sleepWide[i,]~c(0:9))
}
lines(colMeans(sleepWide)~c(0:9),col="blue",lwd=4,lty=2)

## ----eval=FALSE----------------------------------------------------------
## curve(dnorm(x,-1.5),-4.5,3.5)
## xseq = seq(0,3.5,length.out=500)
## yseq = dnorm(xseq,-1.5)
## polygon(x=c(xseq,xseq[length(xseq):1],0),
##         y=c(yseq,rep(0,length(xseq)+1)),
##         col=rgb(0.25,0.75,1,alpha=0.5))

## ----height=5,echo=FALSE,results='hide'----------------------------------
curve(dnorm(x,-1.5),-4.5,3.5)
xseq = seq(0,3.5,length.out=500)
yseq = dnorm(xseq,-1.5)
polygon(x=c(xseq,xseq[length(xseq):1],0),
        y=c(yseq,rep(0,length(xseq)+1)),
        col=rgb(0.25,0.75,1,alpha=0.5))

## ----eval=FALSE----------------------------------------------------------
## text()
## segments()
## arrows()
## symbols()

## ----eval=FALSE----------------------------------------------------------
## jpeg()
## pdf()
## png()
## bmp()
## tiff()

## ----eval=FALSE----------------------------------------------------------
## jpeg("foo.jpg",height=800,width=800)#in pixels
## plot( ... )
## dev.off()

## ----eval=FALSE----------------------------------------------------------
## require("ggplot2")
## blowdown$y <-
##   factor(blowdown$y,labels=c("Surv","Died"))
## ggplot(blowdown,aes(factor(y),D))+
##   geom_boxplot()+
##   theme(axis.text=element_text(size=20),
##         axis.title=element_text(size=0,colour="white"))

## ----echo=FALSE,results='hide'-------------------------------------------
require("ggplot2")
blowdown$y <- 
  factor(blowdown$y,labels=c("Surv","Died"))
ggplot(blowdown,aes(factor(y),D))+
  geom_boxplot()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))

## ----eval=FALSE----------------------------------------------------------
## ggplot(blowdown,aes(D,fill=SPP))+
##   geom_histogram()+
##   theme(axis.text=element_text(size=20),
##         axis.title=element_text(size=0,colour="white"))

## ----echo=FALSE,results='hide'-------------------------------------------
try({
ggplot(blowdown,aes(D,fill=SPP))+
  geom_histogram(bins=30)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))},silent=T)

## ----eval=FALSE----------------------------------------------------------
## ggplot(blowdown,aes(D,fill=SPP,colour=SPP))+
##   geom_density(alpha=0.1)+
##   theme(axis.text=element_text(size=20),
##         axis.title=element_text(size=0,colour="white"))

## ----echo=FALSE,results='hide'-------------------------------------------
ggplot(blowdown,aes(D,fill=SPP,colour=SPP))+
  geom_density(alpha=0.1)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))

## ----eval=FALSE----------------------------------------------------------
## ggplot(blowdown,aes(D,fill=SPP,colour=SPP))+
##   geom_density(position="stack")+
##   theme(axis.text=element_text(size=20),
##         axis.title=element_text(size=0,colour="white"))

## ----echo=FALSE,results='hide'-------------------------------------------
ggplot(blowdown,aes(D,fill=SPP,colour=SPP))+
  geom_density(position="stack")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))

## ----eval=FALSE----------------------------------------------------------
## ggplot(airquality,aes(Solar.R,Ozone))+
##   geom_point(size=6,aes(color=Temp))+
##   geom_smooth(method="lm")+
##   theme(axis.text=element_text(size=20),
##         axis.title=element_text(size=0,colour="white"))

## ----echo=FALSE,results='hide'-------------------------------------------
options(warn=-1)
ggplot(airquality,aes(Solar.R,Ozone))+
  geom_point(size=6,aes(color=Temp))+
  geom_smooth(method="lm")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=0,colour="white"))
options(warn=0)

## ----eval=FALSE----------------------------------------------------------
## normMix = function(x,y){
##   mixProbs = c(1/3,1/2,1/6)
##   ret = dnorm(x,-2)*dnorm(y,-2) +
##     dnorm(x,0)*dnorm(y,2) +
##     dnorm(x,2)*dnorm(y,0)
##   return(ret)
## }
## curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="contour",
##         xlab="",ylab="",labcex=1.5,nlevels=20)
## curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="persp",theta=-15,
##         xlab="",ylab="",zlab="")
## curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="rgl",
##         xlab="",ylab="",zlab="",
##         col = rgb(20/256,60/256,120/256,0.5))

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
require("emdbook")
normMix = function(x,y){
  mixProbs = c(1/3,1/2,1/6)
  ret = dnorm(x,-2)*dnorm(y,-2) +
    dnorm(x,0)*dnorm(y,2) +
    dnorm(x,2)*dnorm(y,0)
  return(ret)
}
curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="contour",
        xlab="",ylab="",labcex=1.5,nlevels=20)

## ----echo=FALSE,results='hide'-------------------------------------------
curve3d(normMix(x,y),from=c(-5,-5),to=c(5,5),sys3d="persp",
        xlab="",ylab="",zlab="",theta=-15)

## ----eval=FALSE----------------------------------------------------------
## require("plot3D");require("alr3")
## pollution = sniffer[,c("TankTemp","TankPres")]
## pollution = data.frame(pollution,logy=log(sniffer$Y))
## with(pollution,scatter3D(TankTemp,TankPres,logy,theta=45,
##                          phi=20,xlab="Temp",ylab="Press",
##                          zlab="Pollution",pch=16))

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
require("plot3D");require("alr3")
pollution = sniffer[,c("TankTemp","TankPres")]
pollution = data.frame(pollution,logy=log(sniffer$Y))
with(pollution,scatter3D(TankTemp,TankPres,logy,theta=45,
                         phi=20,xlab="Temp",ylab="Press",
                         zlab="Pollution",pch=16))

## ----eval=FALSE----------------------------------------------------------
## with(pollution,scatter3D(TankTemp,TankPres,logy,theta=45,
##                          phi=20,xlab="Temp",ylab="Press",
##                          zlab="Pollution",pch=16,
##                          type="h"))

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
with(pollution,scatter3D(TankTemp,TankPres,logy,theta=45,
                         phi=20,xlab="Temp",ylab="Press",
                         zlab="Pollution",pch=16,
                         type="h"))

## ----eval=FALSE----------------------------------------------------------
## fit = lm(logy~TankTemp+TankPres,data=pollution)
## xgrid=with(pollution,
##        seq(min(TankTemp),max(TankTemp),length.out=15))
## ygrid=with(pollution,
##        seq(min(TankPres),max(TankPres),length.out=15))
## xygrid = expand.grid(TankTemp=xgrid,TankPres=ygrid)
## logyPred = matrix(predict(fit,xygrid),15,15)
## fitPts = predict(fit)
## with(pollution,
##      scatter3D(TankTemp,TankPres,logy,theta=45,
##                phi=20,xlab="Temp",ylab="Press",
##                zlab="Pollution",pch=16,
##                surf=list(x=xgrid,y=ygrid,z=logyPred,
##                fit=fitPts,facets=NA)))

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
fit = lm(logy~TankTemp+TankPres,data=pollution)
xgrid=with(pollution,
       seq(min(TankTemp),max(TankTemp),length.out=15))
ygrid=with(pollution,
       seq(min(TankPres),max(TankPres),length.out=15))
xygrid = expand.grid(TankTemp=xgrid,TankPres=ygrid)
logyPred = matrix(predict(fit,xygrid),15,15)
fitPts = predict(fit)
with(pollution,
     scatter3D(TankTemp,TankPres,logy,theta=45,
               phi=20,xlab="Temp",ylab="Press",
               zlab="Pollution",pch=16,
               surf=list(x=xgrid,y=ygrid,z=logyPred,
               fit=fitPts,facets=NA)))

## ----eval=FALSE----------------------------------------------------------
## BDTab = table(cut(blowdown$D,seq(5,85,by=10),
##                   include.lowest=TRUE),
##               blowdown$SPP)
## hist3D(z=BDTab,col = rgb(20/256,60/256,120/256),
##        border = "black",shade = 0.4,space = 0.15,
##        xlab="Diameter",ylab="Species",zlab="")
## plotrgl()

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
BDTab = table(cut(blowdown$D,seq(5,85,by=10),
                  include.lowest=TRUE),
              blowdown$SPP)
hist3D(z=BDTab,col = rgb(20/256,60/256,120/256),
       border = "black",shade = 0.4,space = 0.15,
       xlab="Diameter",ylab="Species",zlab="")

## ----cache=TRUE----------------------------------------------------------
RFun = function(n,start=1){
  ret=0
  for(i in start:n){
    ret = ret + i
  }
  return(ret)
}
RFun(100) == 100*101/2

## ----cache=TRUE----------------------------------------------------------
require("compiler")
RFunComp = cmpfun(RFun)
system.time(RFun(1e7))
system.time(RFunComp(1e7))

## ----cache=TRUE----------------------------------------------------------
require("foreach")
require("doParallel")
registerDoParallel(cl=2)
system.time({
  temp = foreach(i=c(5e6,1e7)) %dopar% {
    if(i==5e6){
      RFun(i)
    }else{
      RFun(i,5e6+1)
    }
  }
  print(sum(unlist(temp)) == 1e7*(1e7+1)/2)
})
stopImplicitCluster()

## ----include=FALSE,cache=TRUE--------------------------------------------
Rcpp::sourceCpp('H:Talks/R_Seminar/RSeminarCppFunctions.cpp')

## ------------------------------------------------------------------------
cppFun(100)

## ----cache=TRUE----------------------------------------------------------
system.time(RFun(1e7))
system.time(cppFun(1e7))

## ----eval=FALSE----------------------------------------------------------
## attach(blowdown)

## ----cache=TRUE----------------------------------------------------------
dim(blowdown)
head(blowdown)

## ----results='hide',echo=FALSE-------------------------------------------
rm(list=setdiff(ls(),"blowdown"))
D <- blowdown$D
S <- blowdown$S
y <- blowdown$y
SPP <- blowdown$SPP

## ----cache=TRUE----------------------------------------------------------
y = factor(y,labels=c("Survived","Died"))
( tab = table(SPP,y) )

## ----cache=TRUE----------------------------------------------------------
chisq.test(tab)
prop.table(tab,1)

## ----eval=FALSE----------------------------------------------------------
## image(t(prop.table(tab,1)),xaxt="n",yaxt="n",
##       main="Tree Species and Survival")
## box()
## axis(1,at=0:1,labels=levels(y),cex.axis=1.5)
## axis(2,at=seq(1,0,length.out=length(levels(SPP))),
##      labels=levels(SPP),las=1,cex.axis=1.5)
## text(cbind(rep(0:1,each=length(levels(SPP))),
##           rep(seq(0,1,length.out=length(levels(SPP))),2)),
##      labels=round(c(prop.table(tab,1)),2),cex=1.5)

## ----echo=F,results='hide',cache=TRUE------------------------------------
image(t(prop.table(tab,1)),xaxt="n",yaxt="n",
      main="Tree Species and Survival")
box()
axis(1,at=0:1,labels=levels(y),cex.axis=1.5)
axis(2,at=seq(1,0,length.out=length(levels(SPP))),
     labels=levels(SPP),las=1,cex.axis=1.5)
text(cbind(rep(0:1,each=length(levels(SPP))),
          rep(seq(0,1,length.out=length(levels(SPP))),2)),
     labels=round(c(prop.table(tab,1)),2),cex=1.5)

## ----cache=TRUE----------------------------------------------------------
logD = log(D)
t.test(logD~y,var.equal=FALSE)

## ----eval=FALSE----------------------------------------------------------
## par(mfrow=c(1,2))
## by(logD,y,qqnorm)
## par(mfrow=c(1,1))

## ----results='hide',echo=FALSE,cache=TRUE--------------------------------
par(mfrow=c(1,2))
by(logD,y,qqnorm)
par(mfrow=c(1,1))

## ----cache=TRUE----------------------------------------------------------
kruskal.test(logD,y)

## ----eval=FALSE----------------------------------------------------------
## ind=which(y=="Died")
## plot(density(x=logD[ind]),ylim=c(0,0.9),
##      main="log(Diameter)")
## lines(density(x=logD[-ind]),lty=2)
## legend("topright",legend=c("Died","Survived"),lty=1:2)

## ----echo=F,results='hide',cache=TRUE------------------------------------
ind=which(y=="Died")
plot(density(x=logD[ind]),ylim=c(0,0.9),
     main="log(Diameter)")
lines(density(x=logD[-ind]),lty=2)
legend("topright",legend=c("Died","Survived"),lty=1:2)

## ----cache=TRUE----------------------------------------------------------
boxplot(logD~SPP)

## ----cache=TRUE----------------------------------------------------------
lmMod= lm(logD~SPP)
anova(lmMod)

## ----cache=TRUE----------------------------------------------------------
summary(lmMod)

## ----eval=FALSE----------------------------------------------------------
## qqnorm(resid(lmMod));qqline(resid(lmMod))
## plot(resid(lmMod)~fitted(lmMod));abline(h=0)

## ----echo=F,results='hide',cache=TRUE------------------------------------
qqnorm(resid(lmMod));qqline(resid(lmMod))
# plot(resid(lmMod)~fitted(lmMod));abline(h=0)

## ----echo=F,results='hide',cache=TRUE------------------------------------
# qqnorm(resid(lmMod));qqline(resid(lmMod))
plot(resid(lmMod)~fitted(lmMod));abline(h=0)

## ----eval=FALSE----------------------------------------------------------
## plot(density(x=S),ylim=c(0,2.25),main="Storm Severity")
## lines(density(x=S[-ind]),lty=2)
## legend("topright",legend=c("Died","Survived"),lty=1:2)

## ----echo=F,results='hide',cache=TRUE------------------------------------
plot(density(x=S),ylim=c(0,2.25),main="Storm Severity")
lines(density(x=S[-ind]),lty=2)
legend("topright",legend=c("Died","Survived"),lty=1:2)

## ----cache=TRUE----------------------------------------------------------
( logRegMod = glm(y~logD + SPP + S, family=binomial) )

## ----cache=TRUE----------------------------------------------------------
round(summary(logRegMod)$coef,4)

## ----cache=TRUE----------------------------------------------------------
( probRegMod = glm(y~logD + SPP + S, 
                   family=binomial(link="probit")) )

## ----cache=TRUE----------------------------------------------------------
summary(probRegMod)$coef

## ----eval=FALSE----------------------------------------------------------
## logitFun = function(x,y){
##   eta = logRegMod$coef["(Intercept)"] +
##     logRegMod$coef["SPPBS"]+logRegMod$coef["logD"]*x+
##     logRegMod$coef["S"]*y
##   return(1/(1+exp(-eta)))
## }
## with(blowdown,
##      curve3d(logitFun,from=c(min(logD),min(S)),
##              to=c(max(logD),max(S)),sys3d="persp",
##              theta=-20,xlab="Diameter",ylab="Storm",
##              zlab="Probability of Dying"))

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
logitFun = function(x,y){
  eta = logRegMod$coef["(Intercept)"] +
    logRegMod$coef["SPPBS"]+logRegMod$coef["logD"]*x+
    logRegMod$coef["S"]*y
  return(1/(1+exp(-eta)))
}
with(blowdown,
     curve3d(logitFun,from=c(min(logD),min(S)),
             to=c(max(logD),max(S)),sys3d="persp",
             theta=-20,xlab="Diameter",ylab="Storm",
             zlab="Probability of Dying"))

## ----eval=FALSE----------------------------------------------------------
## yr <- floor(tt <- time(mdeaths))
## plot(ldeaths,ylab="",
##      xy.labels=paste(month.abb[12*(tt-yr)],
##                      yr-1900,sep="'"))

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
yr <- floor(tt <- time(mdeaths))
plot(ldeaths,ylab="",
     xy.labels=paste(month.abb[12*(tt-yr)],
                     yr-1900,sep="'"))

## ----results='hide'------------------------------------------------------
require("MCMCpack")

## ----cache=TRUE----------------------------------------------------------
require("MCMCpack")
fun1= function(y,nsims=1000,stepsAhead=100){
  TT = length(y)
  BB = sum(y[-TT]^2)
  bb = sum(y[-1]*y[-TT])/BB
  Qb = sum((y[-1]-bb*y[-TT])^2)
  
  s2 = rinvgamma(nsims,shape=0.5*(TT-3),scale=0.5*Qb)
  phi = rnorm(nsims,mean=bb,sd=sqrt(s2/BB))
  ypred = matrix(0.0,nsims,stepsAhead)
  ypred[,1]=rnorm(nsims,phi*y[TT],sd=sqrt(s2))
  for(tt in 2:stepsAhead){
    ypred[,tt] = rnorm(nsims,phi*ypred[,tt-1],sd=sqrt(s2))
  }
  return(list(s2=s2,phi=phi,ypred=ypred))
}

## ----cache=TRUE----------------------------------------------------------
set.seed(1)
M=1000
phiTrue = 0.75
s2True = 0.5
TT=500
Y = matrix(NA,M,TT+1)
Y[,1]=rnorm(M,sd=sqrt(s2True/(1-phiTrue^2)))
phiHat = s2Hat = numeric(M)
for(tt in 1+1:TT){
  Y[,tt] = 0.75*Y[,tt-1] + rnorm(M,sd=sqrt(s2True))
}
for(iter in 1:M){
  temp = fun1(Y[iter,],stepsAhead=2)
  s2Hat[iter] = mean(temp$s2)
  phiHat[iter] = mean(temp$phi)
}

## ----eval=FALSE----------------------------------------------------------
## hist(phiHat)
## abline(v=phiTrue,col="red",lwd=2)
## hist(s2Hat)
## abline(v=s2True,col="red",lwd=2)

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
hist(phiHat)
abline(v=phiTrue,col="red",lwd=2)

## ----echo=FALSE,results='hide',cache=TRUE--------------------------------
hist(s2Hat)  
abline(v=s2True,col="red",lwd=2)

## ----eval=FALSE----------------------------------------------------------
## ldeaths1 = ldeaths-mean(ldeaths)
## fit1 = fun1(ldeaths1,nsims=1e5,stepsAhead=12)
## ypreds = colMeans(fit1$ypred)+mean(ldeaths)
## yBounds = t(apply(fit1$ypred,2,quantile,
##                 probs=c(0.025,0.975)))+mean(ldeaths)
## YL = range(c(ldeaths,yBounds))
## plot(c(ldeaths),ylab="",xlab="",type="l",
##      xlim=c(1,(length(ldeaths)+12)),ylim=YL,lwd=2)
## lines(ypreds~c(length(ldeaths)+1:12),col="red",lwd=2)
## lines(yBounds[,1]~c(length(ldeaths)+1:12),col="blue",
##       lwd=2,lty=2)
## lines(yBounds[,2]~c(length(ldeaths)+1:12),col="blue",
##       lwd=2,lty=2)

## ----echo=FALSE,cache=TRUE-----------------------------------------------
ldeaths1 = ldeaths-mean(ldeaths)
fit1 = fun1(ldeaths1,nsims=1e5,stepsAhead=12)
ypreds = colMeans(fit1$ypred)+mean(ldeaths)
yBounds = t(apply(fit1$ypred,2,quantile,
                probs=c(0.025,0.975)))+mean(ldeaths)
YL = range(c(ldeaths,yBounds))
plot(c(ldeaths),ylab="",xlab="",type="l",
     xlim=c(1,(length(ldeaths)+12)),ylim=YL,lwd=2)
lines(ypreds~c(length(ldeaths)+1:12),col="red",lwd=2)
lines(yBounds[,1]~c(length(ldeaths)+1:12),col="blue",
      lwd=2,lty=2)
lines(yBounds[,2]~c(length(ldeaths)+1:12),col="blue",
      lwd=2,lty=2)

## ----eval=FALSE----------------------------------------------------------
## TT = length(ldeaths1);  BB = sum(ldeaths1[-TT]^2);
## bb = sum(ldeaths1[-1]*ldeaths1[-TT])/BB;
## Qb = sum((ldeaths1[-1]-bb*ldeaths1[-TT])^2)
## curve(dgamma(x,shape=0.5*(TT-3),scale=0.5*Qb),
##       lwd=2,xlab=expression(sigma^2),from=0,to=5e8,
##       ylab="",cex.axis=1.5,cex.lab=1.5)
## curve(dnorm(x,bb,sd=sqrt(mean(fit1$s2)/BB)),
##       lwd=2,xlab=expression(phi),from=0,to=1,
##       ylab="",cex.axis=1.5,cex.lab=1.5)

## ----echo=FALSE,cache=TRUE-----------------------------------------------
TT = length(ldeaths1);  BB = sum(ldeaths1[-TT]^2);  
bb = sum(ldeaths1[-1]*ldeaths1[-TT])/BB;  
Qb = sum((ldeaths1[-1]-bb*ldeaths1[-TT])^2)
curve(dgamma(x,shape=0.5*(TT-3),scale=0.5*Qb),
      lwd=2,xlab=expression(sigma^2),from=0,to=5e8,
      ylab="",cex.axis=1.5,cex.lab=1.5)

## ----echo=FALSE,cache=TRUE-----------------------------------------------
curve(dnorm(x,bb,sd=sqrt(mean(fit1$s2)/BB)),
      lwd=2,xlab=expression(phi),from=0,to=1,
      ylab="",cex.axis=1.5,cex.lab=1.5)

## ----cache=TRUE----------------------------------------------------------
require("compiler")
fun1Comp = cmpfun(fun1)

## ----include=FALSE,cache=TRUE--------------------------------------------
Rcpp::sourceCpp('H:Talks/R_Seminar/RSeminarCppFunctions.cpp')

## ----cache=TRUE----------------------------------------------------------
system.time({for(it in 1:100){
  Rsims = fun1(Y[1,])
  }
})
system.time({for(it in 1:100){
  RsimsComp = fun1Comp(Y[1,])
  }
})
system.time({for(it in 1:100){
  Rcppsims = cppFun1(Y[1,],1000,100)
  }
})

