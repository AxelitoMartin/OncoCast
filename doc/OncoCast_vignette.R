## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.width=9,fig.height=6,message=F,warning = F)
load("../data/survData.rda")
library(survival)
# library(OncoCast)

## ----eval=FALSE----------------------------------------------------------
#  ## SETUP ###
#  set.seed(3)
#  ImpCovMin <- 1
#  ImpCovMax <- 1.5
#  n <- 300
#  
#  x.imp <- cbind(rbinom(n, 1, runif(1,0.1,0.5)),rbinom(n, 1, runif(1,0.1,0.5)),
#                 rbinom(n, 1, runif(1,0.1,0.5)),rbinom(n, 1, runif(1,0.1,0.5)),
#                 rbinom(n, 1, runif(1,0.1,0.5)))
#  
#  coefs <- c(runif(2,ImpCovMin,ImpCovMax),runif(3,-ImpCovMax,-ImpCovMin))
#  x <- x.imp
#  mu <- as.vector(coefs %*% t(x))
#  
#  time <- rexp(n,exp(mu))*12
#  c <- runif(n,0,as.numeric(quantile(time,0.8)))
#  status <- ifelse(time < c , 1,0)
#  time <-pmin(time,c)
#  
#  numDummy <- 10
#  DummyCov <- lapply(1:numDummy,function(x){
#    rbinom(n, 1, runif(1,0.1,0.2))
#  })
#  x.dumb <- do.call(cbind, DummyCov)
#  x <- as.data.frame(cbind(x,x.dumb))
#  
#  data <- as.data.frame(cbind(time,status,x))
#  colnames(data) <- c("time","status",paste0("ImpCov",1:ncol(x.imp)),paste0("Cov",(ncol(x.imp)+1):ncol(x)))
#  rownames(data) <- paste0("Patient",1:n)
#  survData <- data
#  
#  # devtools::use_data(survData, internal = F,overwrite = T)

## ----echo=T--------------------------------------------------------------
### dataset 1
dim(survData)
head(survData[,1:7])

## ------------------------------------------------------------------------
library(OncoCast)

## ------------------------------------------------------------------------
anyNA(survData)
dim(survData)

## ----message=F-----------------------------------------------------------
library(doParallel)
detectCores()

## ------------------------------------------------------------------------
# library(OncoCast)
# library(survival)
Out <- OncoCast(data=survData,formula = Surv(time,status)~.,
                method = "LASSO",runs = 50,
                save = F,nonPenCol = NULL,cores =1)

## ------------------------------------------------------------------------
length(Out$LASSO)
str(Out$LASSO[[1]])

## ----message=F-----------------------------------------------------------
lasso.results <- getResults_OC(Out$LASSO,data=survData,numGroups=4,cuts = c(0.25,0.5,0.75),mut.data = T)

## ----include=F-----------------------------------------------------------
library(knitr)

## ------------------------------------------------------------------------
lasso.results$RiskHistogram
lasso.results$RiskScoreSummary
kable(summary(lasso.results$RiskRefit)$coefficients)
lasso.results$ciSummary

## ----fig.width=9---------------------------------------------------------
lasso.results$inflPlot
lasso.results$selectInflPlot

## ----fig.width=9---------------------------------------------------------
lasso.results$KM

## ----fig.width=9---------------------------------------------------------
lasso.results$mut_Plot

## ----fig.height=7,fig.width=10-------------------------------------------
lasso.results$PieChart

## ------------------------------------------------------------------------
set.seed(3)
new.data <- as.data.frame(matrix(rbinom(10*100,1,0.5),nrow=100,ncol = 10))
colnames(new.data) <- c("ImpCov1","ImpCov2","ImpCov3","ImpCov4","ImpCov5","Cov6","Cov7",
                        "Cov8","Cov9","Cov10")
rownames(new.data) <- paste0("Incoming",1:100)
Incoming <- predIncomingSurv(Out$LASSO,new.data,surv.print = c(5,10,15),riskRefit = lasso.results$RiskRefit)

## ------------------------------------------------------------------------
head(Incoming$data.out)

## ----fig.width=9---------------------------------------------------------
Incoming$RiskHist

## ----fig.width=9---------------------------------------------------------
Incoming$IncKM

## ------------------------------------------------------------------------
# add time component
coefs <- c(1.2,1.3,-1.5,-1,-1,rep(0,ncol(new.data)-5))
mu <- as.vector(coefs %*% t(new.data))
n = nrow(row(new.data))
time <- rexp(n,exp(mu))*12
c <- runif(n,0,as.numeric(quantile(time,0.8)))
status <- ifelse(time < c , 1,0)
time <-pmin(time,c)
new.data$time <- time
new.data$status <- status

## ------------------------------------------------------------------------
validation <- validate(OC_object=Out$LASSO,Results=lasso.results,
in.data=new.data, formula=Surv(time,status)~.)

## ----fig.width=9---------------------------------------------------------
validation$RiskHistogram.new

## ------------------------------------------------------------------------
validation$out.data[1:10,c("OncoCastRiskScore","RiskGroup")]

## ----fig.width=9---------------------------------------------------------
validation$KM

## ------------------------------------------------------------------------
head(survData.LT[,1:5]) 

## ------------------------------------------------------------------------
# launchApp()

## ----pressure, echo=FALSE, out.width = '100%'----------------------------
knitr::include_graphics("preview.png")

## ---- echo=FALSE, out.width = '100%'-------------------------------------
knitr::include_graphics("results1.png")

## ---- echo=FALSE, out.width = '100%'-------------------------------------
knitr::include_graphics("results2.png")

## ---- echo=FALSE, out.width = '100%'-------------------------------------
knitr::include_graphics("results3.png")

