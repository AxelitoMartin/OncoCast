#' validate
#'
#' This function let's user validate an OncoCast model generated in a new cohort which has it's own
#' survival data.
#'
#' @param OC_object Output of the OncoCast function.
#' @param Fits All fits from the OncoCast run.
#' @param new.data New data set containing the information of the incoming patients. Should be a dataframe
#' with patients as rows and features as columns.
#' @param limit Optional numerical argument to set a time limit on the KM plot
#' @return data.out : The data frame inputted in the function with an additional column giving the predicted risk
#' score of the incoming patients.
#' @return RiskHist : A histogram of the distribution of the risk scores of patients in the given dataset.
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula = Surv(time,status)~.,
#' method = "LASSO",runs = 50,
#' save = F,nonPenCol = NULL,cores =2)
#' results <- getResults_OC(OC_object=test$LASSO,data=survData,
#' numGroups=5,cuts=c(0.2,0.4,0.6,0.8),
#' geneList="NULL",mut.data=T)
#' new.data <- as.data.frame(matrix(rbinom(5*20,1,0.5),nrow=20,ncol = 5))
#' colnames(new.data) <- c("ImpCov1","ImpCov2","ImpCov3","ImpCov4","Cov7")
#' rownames(new.data) <- paste0("Incoming",1:20)
#'
#' # add time component
#' coefs <- c(1.2,1.3,-1.5,-1,-1,rep(0,ncol(new.data)-5))
#' mu <- as.vector(coefs %*% t(new.data))
#' n = nrow(row(new.data))
#' time <- rexp(n,exp(mu))*12
#' c <- runif(n,0,as.numeric(quantile(time,0.8)))
#' status <- ifelse(time < c , 1,0)
#' time <-pmin(time,c)
#' new.data$time <- time
#' new.data$status <- status
#'
#' validation <- validate(OC_object=test$LASSO,Results=results,
#' in.data=new.data, formula=Surv(time,status)~.)



validate <- function(OC_object,Results,in.data,formula,limit = NULL){

  OC_object <- Filter(Negate(is.null), OC_object)
  means.train <- sapply(OC_object,"[[","means")

  qts = Results$rawCuts
  ori.risk <- as.numeric(Results$risk.raw)

  if(OC_object[[1]]$method %in% c("LASSO","RIDGE","ENET")){

    LassoFits <- as.matrix(Results$Fits)
    features <- colnames(LassoFits)

    if(!all(is.na(match(colnames(in.data),features)))){
      matched.genes <- c(na.omit(match(colnames(in.data),features)))
      new.dat <- in.data[,which(!is.na(match(colnames(in.data),features)))]

      ## ADD ALL MISSING GENES TO BE ALL zero ##
      missing <- features[which(is.na(match(features,colnames(new.dat))))]
      to.add <- as.data.frame(matrix(0L,nrow=nrow(new.dat),ncol=length(missing)))
      colnames(to.add) <- missing
      rownames(to.add) <- rownames(new.dat)
      new.dat <- as.data.frame(cbind(new.dat,to.add))

      new.dat <- new.dat[,match(features,colnames(new.dat))]

      #############################################

      all.pred <- lapply(1:nrow(LassoFits),function(x){

        ### Subset to the coefs of that cv ###
        coefs <- LassoFits[x,LassoFits[x,] != 0]
        new.temp <- select(new.dat,names(coefs))

        ## substract mean mutation rate of TRAINING SET !!!###
        new.x <- new.temp - rep(means.train[[x]][match(names(coefs),names(means.train[[x]]))], each = nrow(new.temp))
        cal.risk.test <- drop(as.matrix(new.x) %*% coefs)
        return(cal.risk.test)
      })

    }
    else{
      stop("No gene in your dataset overlapped with the original data. Please rename genes or check your dataset.")
    }

  }


  if(OC_object[[1]]$method %in% c("GBM","RF","SVM")){

    if(OC_object[[1]]$method == "GBM") features <- OC_object[[1]]$GBM$var.names
    if(OC_object[[1]]$method == "RF") features <- rownames(OC_object[[1]]$RF$importance)
    if(OC_object[[1]]$method == "SVM") features <- names(OC_object[[1]]$Vars)
    if(!all(is.na(match(colnames(in.data),features)))){
      matched.genes <- c(na.omit(match(colnames(in.data),features)))
      new.dat <- in.data[,which(!is.na(match(colnames(in.data),features)))]

      ## ADD ALL MISSING GENES TO BE ALL zero ##
      missing <- features[which(is.na(match(features,colnames(new.dat))))]
      to.add <- as.data.frame(matrix(0L,nrow=nrow(new.dat),ncol=length(missing)))
      colnames(to.add) <- missing
      rownames(to.add) <- rownames(new.dat)
      new.dat <- as.data.frame(cbind(new.dat,to.add))

      new.dat <- new.dat[,match(features,colnames(new.dat))]

      if(OC_object[[1]]$method == "GBM") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$GBM,newdata=new.dat,
                  n.trees = x$bestTreeForPrediction,
                  type="response")
        })}

      if(OC_object[[1]]$method == "RF") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$RF,new.dat)
        })}

      if(OC_object[[1]]$method == "SVM") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$SVM,new.dat)
        })}
    }
    else{
      stop("No gene in your dataset overlapped with the original data. Please rename genes or check your dataset.")
    }
  }


  ### PUTTING THINGS TOGETHER ####
  all.pred <- do.call("cbind",all.pred)
  Risk <- apply(all.pred,1,mean)
  names(Risk) <- rownames(new.dat)
  # Risk.all <- as.matrix(coefs) %*% as.matrix(t(new.dat))
  # Risk <- apply(Risk.all,2,mean)
  #in.data$Risk <- Risk
  ##########################################
  ori.risk.range <- range(ori.risk)
  in.data$OncoCastRiskScore <- rescale(Risk, to = c(0, 10), from = ori.risk.range) #WithOriginal
  #in.data$rescaledRisk <- rescale(in.data$Risk, to = c(0, 10), from = range(in.data$Risk, na.rm = TRUE, finite = TRUE))
  RiskHistogram.new <- ggplot(in.data, aes(x = OncoCastRiskScore, y = ..density..)) +
    geom_histogram(show.legend = FALSE, aes(fill=..x..),
                   breaks=seq(min(in.data$OncoCastRiskScore,na.rm = T), max(in.data$OncoCastRiskScore,na.rm = T), by=20/nrow(in.data))) +
    geom_density(show.legend = FALSE) +
    theme_minimal() +
    labs(x = "Average risk score", y = "Density") +
    scale_fill_gradient(high = "red", low = "green")

  #return(list("RiskHistogram.new"=RiskHistogram.new,"out.data"=in.data))


  #### NEED TO REMAKE CUTS BASED ON THE ORIGINAL SCORE ####

  #### deal with surv formula #####
  survFormula <- as.formula(formula)
  survResponse <- survFormula[[2]]

  ### reprocess data
  if(length(as.list(survResponse)) == 3){
    colnames(in.data)[match(as.list(survResponse)[2:3],colnames(in.data))] <- c("time","status")
    LT = FALSE
  }
  if(length(as.list(survResponse)) == 4){
    colnames(in.data)[match(as.list(survResponse)[2:4],colnames(in.data))] <- c("time1","time2","status")
    LT = TRUE
  }

  riskGroup <- c()
  numGroups <- length(qts) + 1
  for(i in 1:length(Risk)){

    if(numGroups == 2){
      if(Risk[i] > qts) riskGroup[i] <- 2
      else riskGroup[i] <- 1
    }

    if(numGroups >2){
      temp <- Risk[i] - qts
      ind1 <- which(temp>0)
      ind2 <- which(temp<0)
      if(length(ind1)>0 && length(ind2)>0) {riskGroup[i] <- ind2[1]}
      else if (length(ind1)==0){riskGroup[i] <- 1}
      else if (length(ind2)==0){riskGroup[i] <- numGroups}
    }
  }

  # if(numGroups == 2){riskGroup[is.na(riskGroup)] <- 1}
  in.data$RiskGroup <- riskGroup
  in.data$RiskGroup <- factor(in.data$RiskGroup, levels = c(1:numGroups) )


  if(LT == TRUE) {
    fit0 <- coxph(Surv(time1,time2,status) ~ RiskGroup,data=in.data,
                  na.action=na.exclude)
    if(max(in.data$time2) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  if(LT == FALSE) {
    fit0 <- coxph(Surv(time,status) ~ RiskGroup,data=in.data,
                  na.action=na.exclude)
    if(max(in.data$time) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }

  log.test.pval <- as.vector(summary(fit0)[10][[1]])[3]
  CI <- as.numeric(as.vector(summary(fit0)[14])[[1]][1])

  if(is.null(limit)) {
    if(LT) limit <- as.numeric(max(in.data$time2))
    if(!LT) limit <- as.numeric(max(in.data$time))
  }

  if(length(qts) == 2){
    if(LT) {KM <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,#surv.median.line = "hv",
                             data = in.data,break.time.by = 6,xlim=c(0,limit),risk.table = T,palette = c("chartreuse3","cyan4","darkgoldenrod2")) + xlab("Time") +
      labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))} #,xlim=c(0,limit)
    if(!LT){KM <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,#surv.median.line = "hv",
                             data = in.data,break.time.by = 6,xlim=c(0,limit),risk.table = T,palette = c("chartreuse3","cyan4","darkgoldenrod2")) + xlab("Time") +
      labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))}
  }
  else{
    if(LT) {KM <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,#surv.median.line = "hv",
                             data = in.data,break.time.by = 6,xlim=c(0,limit),risk.table = T) + xlab("Time") +
      labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))} #,xlim=c(0,limit)
    if(!LT){KM <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,#surv.median.line = "hv",
                             data = in.data,break.time.by = 6,xlim=c(0,limit),risk.table = T) + xlab("Time") +
      labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))}
  }

  return(list("RiskHistogram.new"=RiskHistogram.new,"out.data"=in.data,"KM"=KM))

}

