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



validate <- function(OC_object,Results,in.data,formula,limit = NULL,...){

  args <- list(...)
  surv.median.line <- ifelse(is.null(args[['surv.median.line']]),"hv",args[['surv.median.line']])
  risk.table <- ifelse(is.null(args[['risk.table']]),T,args[['risk.table']])

  OC_object <- Filter(Negate(is.null), OC_object)
  means.train <- sapply(OC_object,"[[","means")
  qts = Results$rawCuts
  ori.risk <- as.numeric(Results$risk.raw)

  if(OC_object[[1]]$method %in% c("LASSO","RIDGE","ENET")){

    LassoFits <- as.matrix(Results$Fits)
    features <- colnames(LassoFits)

    dums <- apply(in.data,2,function(x){anyNA(as.numeric(as.character(x)))})
    if(sum(dums) > 0){
      tmp <- in.data %>%
        select(which(dums)) %>%
        fastDummies::dummy_cols(remove_first_dummy = T) %>%
        select(-one_of(names(which(dums))))
      in.data <- as.data.frame(cbind(
        in.data %>% select(-one_of(names(which(dums)))),
        tmp
      ) %>% mutate_all(as.character) %>%
        mutate_all(as.numeric)
      )
      warning("Character variables were transformed to dummy numeric variables. If you didn't have any character variables make sure all columns in your input data are numeric. The transformed data will be saved as part of the output.")
    }

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


  if(OC_object[[1]]$method %in% c("GBM","RF","SVM","NN")){

    if(OC_object[[1]]$method == "GBM") features <- OC_object[[1]]$GBM$var.names
    if(OC_object[[1]]$method == "RF") features <- OC_object[[1]]$RF$forest$independent.variable.names
    if(OC_object[[1]]$method == "SVM") features <- names(OC_object[[1]]$Vars)
    if(OC_object[[1]]$method == "NN") features <- names(OC_object[[1]]$Vars)

    dums <- apply(in.data,2,function(x){anyNA(as.numeric(as.character(x)))})
    if(sum(dums) > 0){
      tmp <- in.data %>%
        select(which(dums)) %>%
        fastDummies::dummy_cols(remove_first_dummy = T) %>%
        select(-one_of(names(which(dums))))
      in.data <- as.data.frame(cbind(
        in.data %>% select(-one_of(names(which(dums)))),
        tmp
      ) %>% mutate_all(as.character) %>%
        mutate_all(as.numeric)
      )
      warning("Character variables were transformed to dummy numeric variables. If you didn't have any character variables make sure all columns in your input data are numeric. The transformed data will be saved as part of the output.")
    }

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
          predict(x$RF,new.dat)$predictions
        })}

      if(OC_object[[1]]$method == "SVM") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$SVM,new.dat)
        })}

      if(OC_object[[1]]$method == "NN") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$NN,new.dat)
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

  if(LT) {KM <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = surv.median.line,
                           data = in.data,break.time.by = 6,xlim=c(0,limit),risk.table = risk.table) + xlab("Time") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))} #,xlim=c(0,limit)
  if(!LT){KM <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = surv.median.line,
                           data = in.data,break.time.by = 6,xlim=c(0,limit),risk.table = risk.table) + xlab("Time") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))}


  survivalGroup <- as.data.frame(matrix(nrow=numGroups,ncol=4))
  rownames(survivalGroup) <- paste0("riskGroup ",1:numGroups)
  colnames(survivalGroup) <- c("MedianOS","95%CI","1Ysurvival","3Ysurvival")
  # for each group find closest value to median
  if(timeType == "Months"){YR1 <- 1*12;YR3 <- 3*12}
  if(timeType == "Days"){YR1 <- 1*365;YR3 <- 3*365}
  for(i in 1:numGroups){
    if(LT == TRUE){NewObject <- with(in.data[in.data$RiskGroup == i,],Surv(time1,time2,status))}
    if(LT == FALSE){NewObject <- with(in.data[in.data$RiskGroup == i,],Surv(time,status))}
    Fit <- survfit(NewObject ~ 1,data=in.data[in.data$RiskGroup == i,], conf.type = "log-log")
    # med.index <- which.min(abs(Fit$surv-0.5))
    YR1.index <- which.min(abs(Fit$time-YR1))
    YR3.index <- which.min(abs(Fit$time-YR3))
    survivalGroup[i,] <- c(as.numeric(round(summary(Fit)$table[7],digits=2)),
                           paste0("(",as.numeric(round(summary(Fit)$table[8],digits=2)),",",
                                  as.numeric(round(summary(Fit)$table[9],digits=2)),")"),
                           paste0(round(Fit$surv[YR1.index],digits=2)," (",
                                  round(Fit$lower[YR1.index],digits=2),",",
                                  round(Fit$upper[YR1.index],digits=2),")"),
                           paste0(round(Fit$surv[YR3.index],digits=2)," (",
                                  round(Fit$lower[YR3.index],digits=2),",",
                                  round(Fit$upper[YR3.index],digits=2),")"))
  }

  return(list("RiskHistogram.new"=RiskHistogram.new,"out.data"=in.data,"KM"=KM,"survTable"=survivalGroup))

}

