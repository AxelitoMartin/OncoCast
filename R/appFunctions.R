#######################################
##### EXTRA FUNCTIONS FOR THE APP #####
#######################################

#### LOAD FUNCTION ####
#' load_object
#' @export
load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

##### get CI and other basic stuff #####
#' Progno.tab
#' @export
Progno.tab <- function(OC_object,data,geneList=NULL){

  method <- OC_object[[1]]$method

  ## determine if left truncated
  if(length(grep("time",colnames(data)))  == 1) {LT = FALSE}
  if(length(grep("time",colnames(data)))  == 2) {LT = TRUE}

  MD <- 12

  ConcordanceIndex <- as.data.frame(as.vector(unlist(sapply(OC_object, "[[", "CI"))))
  summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0.1,0.25,0.5,0.75,0.9),na.rm = T))),digits = 2)
  colnames(summary.CI) <- "Concordance Index"
  rownames(summary.CI) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
  CI.BP <- as.data.frame(t(summary.CI))


  final.pred <- sapply(OC_object,"[[","predicted")
  average.risk <- apply(final.pred,1,function(x){
    mean(as.numeric(x),na.rm = TRUE)
  })
  average.risk[which(is.na(average.risk))] <- NA
  to <- c(0,10)
  from <- range(average.risk, na.rm = TRUE, finite = TRUE)
  RiskScore <- (as.numeric(average.risk)-from[1])/diff(from)*diff(to)+to[1]
  #RiskScore <- rescale(as.numeric(average.risk), to = c(0, 10), from = range(average.risk, na.rm = TRUE, finite = TRUE))
  summary.RiskScore <- round(as.data.frame(c(quantile(RiskScore,c(0.1,0.25,0.33,0.5,0.66,0.75,0.9),na.rm = TRUE))),digits = 2)
  colnames(summary.RiskScore) <- "Risk Score"
  rownames(summary.RiskScore) <- c("Lower 10%","1st Quarter","1st Tertile","Median","2nd Tertile","3rd Quarter","Upper 10%")
  ## refit coxph model with average risk as covariate
  meanRS <- mean(RiskScore)
  #RiskScore <- average.risk #- meanRS
  if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~RiskScore)
  if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~RiskScore)
  Risk <- as.data.frame(RiskScore)

  RiskHistogram <- ggplot(Risk, aes(x = RiskScore, y = ..density..)) +
    geom_histogram(show.legend = FALSE, aes(fill=..x..),
                   breaks=seq(min(Risk$RiskScore,na.rm = T), max(Risk$RiskScore,na.rm = T), by=0.25)) +
    geom_density(show.legend = FALSE) +
    theme_minimal() +
    labs(x = "Average risk score", y = "Density") +
    scale_fill_gradient(high = "red", low = "green")


  return(list("CI" = CI.BP,"risk.raw"=average.risk,"scaled.risk"=RiskScore,
              "RiskHistogram"=RiskHistogram,"RiskScoreSummary"=as.data.frame(t(summary.RiskScore)),
              "RiskRefit"=refit.risk,"LT"=LT))
}

#######################################################


#######################################################

##### RISK STRAT FUNCTION #####
#' riskStrat
#' @export

riskStrat <- function(OC_object,data,cuts,average.risk,LT,plotQuant=1){

  cuts <- as.numeric(cuts)/100
  if(max(cuts) >= 1 || min(cuts) <= 0){
    stop("Inappropriate bounds for stratification")
  }

  ######### RISK STRATIFICATION ############
  # generate groups #
  #if(length(cuts) != (numGroups-1)){stop("Mismatch between number of groups and cuts! Length of cuts argument should equal numGroups-1")}
  numGroups = length(cuts) +1
  qts <- quantile(average.risk,cuts)
  riskGroup <- c()
  for(i in 1:length(average.risk)){
    temp <- average.risk[i] - qts
    if(length(match(names(which.max(temp[temp<0])),names(qts))) == 1) riskGroup[i] <- match(names(which.max(temp[temp<0])),names(qts))
    else riskGroup[i] <- numGroups
  }

  if(numGroups == 2){riskGroup[is.na(riskGroup)] <- 1}
  data$RiskGroup <- riskGroup
  data$RiskGroup <- factor(data$RiskGroup, levels = c(1:numGroups) )

  if(LT == TRUE) {
    fit0 <- coxph(Surv(time1,time2,status) ~ RiskGroup,data=data,
                  na.action=na.exclude)
    if(max(data$time2) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  if(LT == FALSE) {
    fit0 <- coxph(Surv(time,status) ~ RiskGroup,data=data,
                  na.action=na.exclude)
    if(max(data$time) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }

  log.test.pval <- as.vector(summary(fit0)[10][[1]])[3]
  CI <- as.numeric(as.vector(summary(fit0)[14])[[1]][1])
  if(LT) limit <- as.numeric(quantile(data$time2,plotQuant))
  if(!LT) limit <- as.numeric(quantile(data$time,plotQuant))

  if(LT) {KM <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                           data = data,xlim=c(0,limit),break.time.by = 6,risk.table=T) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =4),")",sep=""))}
  if(!LT){KM <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                           data = data,xlim=c(0,limit),break.time.by = 6,risk.table=T) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =4),")",sep=""))}



  survivalGroup <- as.data.frame(matrix(nrow=numGroups,ncol=4))
  rownames(survivalGroup) <- 1:numGroups
  colnames(survivalGroup) <- c("MedianOS","95%CI","1Ysurvival","3Ysurvival")
  # for each group find closest value to median
  if(timeType == "Months"){YR1 <- 1*12;YR3 <- 3*12}
  if(timeType == "Days"){YR1 <- 1*365;YR3 <- 3*365}
  for(i in 1:numGroups){
    if(LT == TRUE){NewObject <- with(data[data$RiskGroup == i,],Surv(time1,time2,status))}
    if(LT == FALSE){NewObject <- with(data[data$RiskGroup == i,],Surv(time,status))}
    Fit <- survfit(NewObject ~ 1,data=data[data$RiskGroup == i,], conf.type = "log-log")
    # med.index <- which.min(abs(Fit$surv-0.5))
    YR3.index <- which.min(abs(Fit$time-YR1))
    YR5.index <- which.min(abs(Fit$time-YR3))
    survivalGroup[i,] <- c(as.numeric(round(summary(Fit)$table[7],digits=2)),
                           paste0("(",as.numeric(round(summary(Fit)$table[8],digits=2)),",",
                                  as.numeric(round(summary(Fit)$table[9],digits=2)),")"),
                           round(Fit$surv[YR3.index],digits=2),round(Fit$surv[YR5.index],digits=2))
  }


  return(list("RiskGroup"=riskGroup,"KM"=KM,"survivalTable"=survivalGroup,"rawCuts"= as.numeric(qts)))

}

########################################


########################################
#' gene.view
#' @export

gene.view <- function(OC_object,data,geneList=NULL,LT,average.risk){

  method <- OC_object[[1]]$method
  ########################
  if(method %in% c("LASSO","RIDGE","ENET")){
    allCoefs <- t(sapply(OC_object,"[[","fit"))
    allCoefs[is.na(allCoefs)] <- 0

    #selected.genes.lasso <- apply(allCoefs,2,function(x){sum(x!=0)})

    #topHits <- names(sort(selected.genes.lasso[selected.genes.lasso >0.5],decreasing = TRUE))
    # if(length(selected.genes.lasso)  >= 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:20]}
    # if(length(selected.genes.lasso)  < 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:length(selected.genes.lasso)]}
    # melt.rank.lasso <- melt(ranked.lasso)
    # melt.rank.lasso$Gene <- factor(rownames(melt.rank.lasso), levels = rownames(melt.rank.lasso))
    # melt.rank.lasso$value <- melt.rank.lasso$value/length(OC_object)
    # colnames(melt.rank.lasso) <- c("Frequency","Gene")
    #
    # influencePlot <- ggplot(melt.rank.lasso,aes(x=Gene,y=Frequency,fill=Gene))+geom_col()+
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #   labs(title = paste("Most recurrent selected genes out of",length(OC_object),"runs")) +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")

    # if(LT) Variables <- colnames(data)[-c(1:3)]
    # if(!LT) Variables <- colnames(data)[-c(1:2)]

    meanCoefs <- apply(allCoefs,2,function(x){mean(x,na.rm = TRUE)})
    selectFreq <- apply(allCoefs,2,function(x){
      length(which(x!=0))/length(x)
    })

    if(length(selectFreq[selectFreq > 0.5]) > 2) {
      topHits <- names(selectFreq[order(selectFreq,decreasing = T)])[1:sum(selectFreq > 0.5)] }

    else topHits <- names(selectFreq[order(selectFreq,decreasing = T)])[1:10]
    ## get mu freq
    if(LT) data.temp <- data[,-c(1:3)]
    if(!LT) data.temp <- data[,-c(1:2)]
    MutationFrequency <- apply(data.temp,2,function(x){
      sum(x)/length(x)
    })

    resultsAll <- as.data.frame(cbind(meanCoefs,selectFreq,MutationFrequency))
    colnames(resultsAll) <- c("MeanCoefficient","SelectionFrequency","MutationFrequency")
    rownames(resultsAll) <- names(meanCoefs)
    resultsAll <- resultsAll[complete.cases(resultsAll),]
    resultsAll$GeneName <- rownames(resultsAll)
    resultsAll$MutationFrequency2 <- cut(resultsAll$MutationFrequency, c(0,0.10,0.20,0.40))

    if(length(geneList)!=0){
      m <- resultsAll[match(geneList,rownames(resultsAll)), ]

      a <- list(
        x = m$MeanCoefficient,
        y = m$SelectionFrequency,
        text = rownames(m),
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 7,
        ax = 20,
        ay = -40
      )

      selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                                text = ~paste('Gene :',GeneName,
                                              '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                                mode = "markers") %>% #,size = ~MutationFrequency,color = ~MutationFrequency
        layout(title ="Volcano Plot",annotations = a)
    }

    else{
      selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                                text = ~paste('Gene :',GeneName,
                                              '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                                mode = "markers") %>% #,size = ~MutationFrequency,color = ~MutationFrequency) %>%
        layout(title ="Volcano Plot")
    }

  }


  if(method %in% c("GBM","RF","SVM")) {
    #### Variables ####
    # if(LT) Variables <- colnames(data)[-c(1:3)]
    # if(!LT) Variables <- colnames(data)[-c(1:2)]

    imp <- sapply(OC_object,"[[","Vars")
    #imp <- do.call("cbind",OC_object$Vars)
    mean.imp <- apply(imp,1,mean)

    # for 10 top medians make boxplot
    topHits <- sort(mean.imp,decreasing = T)[1:pmin(ncol(data)-2,20)]
    topHitsMat <- t(imp[match(names(topHits),rownames(imp)),])
    topHitsMat.melt <- melt(topHitsMat)[,-1]
    colnames(topHitsMat.melt) <- c("Gene","Rank")

    p <- ggplot(topHitsMat.melt, aes(x=Gene,y = Rank)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
      #geom_jitter(width = 0.05) +
      theme_grey() +
      labs(title = "Variable use in the forest", subtitle = "Using genetic data only",y = "Importance") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    selectInflPlot <- ggplotly(p)

    allCoefs <- NULL
  }


  ##### heatmap #####
  average.risk.order <- order(average.risk,decreasing = F)
  risk.ordered <- average.risk[average.risk.order]

  if(LT) uniques <- apply(data[,-c(1:3)],2,function(x){length(unique(x))})
  else uniques <- apply(data[,-c(1:2)],2,function(x){length(unique(x))})

  ## for binary ##
  heatmap.sorted.bin <- NULL

  # data$Stage3 <- ifelse(data$Stage ==2,1,0)
  # data$Stage2 <- ifelse(data$Stage ==1,1,0)

  if(sum(uniques == 2) > 2){
    genes <- names(uniques[which(uniques == 2)])
    if(method %in% c("LASSO","RIDGE","ENET") ) keep <- names(sort(selectFreq[match(genes,names(selectFreq))],decreasing = T)[1:pmin(15,length(genes))])
    if(method %in% c("GBM","RF","SVM")) keep <- names(sort(mean.imp[match(genes,names(mean.imp))],decreasing = T)[1:pmin(15,length(genes))])

    if(!is.null(geneList)){
      if(length(geneList) < 15){
        keep <- c(geneList,keep[-c((15-length(geneList)):15)])
      }
    }

    data.sub <- data[,match(keep,colnames(data))]
    data.ordered <- t(as.matrix(data.sub[average.risk.order,]))
    heatmap.sorted.bin <- pheatmap(data.ordered,cluster_rows = F,cluster_cols = F,show_rownames=T,show_colnames = F,
                                   color = c("cyan","black"),border_color = "black",silent =T,legend=F,
                                   fontsize = 15)
  }

  ## for continuos ##
  heatmap.sorted.cont <- NULL

  if(sum(uniques > 2) > 2){
    genes <- names(uniques[which(uniques > 2)])
    if(method %in% c("LASSO","RIDGE","ENET") ) keep <- names(sort(selectFreq[match(genes,names(selectFreq))],decreasing = T)[1:pmin(15,length(genes))])
    if(method %in% c("GBM","RF","SVM")) keep <- names(sort(mean.imp[match(genes,names(mean.imp))],decreasing = T)[1:pmin(15,length(genes))])

    if(!is.null(geneList)){
      if(length(geneList) < 15){
        keep <- c(geneList,keep[-c((15-length(geneList)):15)])
      }
    }

    data.sub <- data[,match(keep,colnames(data))]
    data.ordered <- apply(data.sub,2,function(x){
      y <- (x-min(x))/max(x-min(x))
    })
    data.ordered.new <- t(as.matrix(data.ordered[average.risk.order,]))
    # data.ordered <- t(as.matrix(data.sub[average.risk.order,]))
    # data.ordered <- scale(data.ordered)
    # data.ordered.new <- t(apply(data.ordered,1,function(x){
    #   to <- c(0,1)
    #   from <- range(x, na.rm = TRUE, finite = TRUE)
    #   new <- (as.numeric(x)-from[1])/diff(from)*diff(to)+to[1]
    #   return(new)
    # }))
    heatmap.sorted.cont <- pheatmap(data.ordered.new,cluster_rows = F,cluster_cols = F,show_rownames=T,show_colnames = F,silent =T)
  }

  return(list("Fits"=allCoefs,"selectInflPlot" = selectInflPlot,
              "heatmap.sorted.bin"=heatmap.sorted.bin,"heatmap.sorted.cont"=heatmap.sorted.cont,"topHits" = topHits))

}


#########################################################


#########################################################
#' validate.app
#' @export

validate.app <- function(OC_object,LassoFits,ori.risk,qts,in.data,formula){

  # deal with formula #
  formula <- gsub(" ","",formula)
  if(formula == "") formula <- OC_object[[1]]$formula
  else{
    parts <- unlist(strsplit(formula,","))
    if(length(parts) == 2) formula <- as.formula(paste0("Surv(",parts[1],",",parts[2],")~."))
    if(length(parts) == 3) formula <- as.formula(paste0("Surv(",parts[1],",",parts[2],",",parts[2],")~."))
  }

  OC_object <- Filter(Negate(is.null), OC_object)
  ori.risk <- as.numeric(ori.risk)

  ### PENALIZED ###

  if(OC_object[[1]]$method %in% c("LASSO","RIDGE","ENET")){
    means.train <- sapply(OC_object,"[[","means")

    LassoFits <- as.matrix(LassoFits)
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

  ### GBM ###

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
    temp <- Risk[i] - qts
    ind1 <- which(temp>0)
    ind2 <- which(temp<0)
    if(length(ind1)>0 && length(ind2)>0) {riskGroup[i] <- ind2[1]}
    else if (length(ind1)==0){riskGroup[i] <- 1}
    else if (length(ind2)==0){riskGroup[i] <- numGroups}
    # if(length(match(names(which.max(temp[temp<0])),names(qts))) == 1) riskGroup[i] <- match(names(which.max(temp[temp<0])),names(qts))
    # else riskGroup[i] <- numGroups
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


  if(LT) {KM <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                           data = in.data,break.time.by = 6) + xlab("Time") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))} #,xlim=c(0,limit)
  if(!LT){KM <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=in.data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                           data = in.data,break.time.by = 6) + xlab("Time") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))}


  return(list("RiskHistogram.new"=RiskHistogram.new,"out.data"=in.data,"KM"=KM))

}


######################################################


######################################################

#' predIncomingSurv.app
#' @export


predIncomingSurv.app <- function(OC_object,new.data,surv.print= NULL,riskRefit){

  OC_object <- Filter(Negate(is.null), OC_object)
  # get all information needed from the oncocast object
  # 1. risk
  final.pred <- sapply(OC_object,"[[","predicted")
  ori.risk <- apply(final.pred,1,function(x){
    mean(as.numeric(x),na.rm = TRUE)
  })

  # 2. Fits


  if(OC_object[[1]]$method %in% c("LASSO","RIDGE","ENET")){
    LassoFits <- t(sapply(OC_object,"[[","fit"))
    LassoFits[is.na(LassoFits)] <- 0

    ################################
    features <- colnames(LassoFits)


    if(!all(is.na(match(colnames(new.data),features)))){
      matched.genes <- c(na.omit(match(colnames(new.data),features)))
      new.dat <- new.data[,which(!is.na(match(colnames(new.data),features)))]

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
        new.x <- new.temp - rep(OC_object[[x]]$means[match(names(coefs),names(OC_object[[x]]$means))], each = nrow(new.temp))
        cal.risk.test <- drop(as.matrix(new.x) %*% coefs)
        return(cal.risk.test)
      })
    }

  }

  if(OC_object[[1]]$method %in% c("GBM","RF","SVM")){

    if(OC_object[[1]]$method == "GBM") features <- OC_object[[1]]$GBM$var.names
    if(OC_object[[1]]$method == "RF") features <- rownames(OC_object[[1]]$RF$importance)
    if(OC_object[[1]]$method == "SVM") features <- names(OC_object[[1]]$Vars)

    if(!all(is.na(match(colnames(new.data),features)))){
      matched.genes <- c(na.omit(match(colnames(new.data),features)))
      new.dat <- new.data[,which(!is.na(match(colnames(new.data),features)))]

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

  all.pred <- do.call("cbind",all.pred)
  Risk <- apply(all.pred,1,mean)
  names(Risk) <- rownames(new.dat)
  # Risk.all <- as.matrix(coefs) %*% as.matrix(t(new.dat))
  # Risk <- apply(Risk.all,2,mean)
  #new.data$Risk <- Risk
  ##########################################
  ori.risk.range <- range(ori.risk)
  new.data$OncoCastRiskScore <- rescale(Risk, to = c(0, 10), from = ori.risk.range) #WithOriginal
  #new.data$rescaledRisk <- rescale(new.data$Risk, to = c(0, 10), from = range(new.data$Risk, na.rm = TRUE, finite = TRUE))
  RiskHistogram.new <- ggplot(new.data, aes(x = OncoCastRiskScore, y = ..density..)) +
    geom_histogram(show.legend = FALSE, aes(fill=..x..),
                   breaks=seq(min(new.data$OncoCastRiskScore,na.rm = T), max(new.data$OncoCastRiskScore,na.rm = T))) +#, by=20/nrow(new.data))) +
    geom_density(show.legend = FALSE) +
    theme_minimal() +
    labs(x = "Average risk score", y = "Density") +
    scale_fill_gradient(high = "red", low = "green")

  #return(list("RiskHistogram.new"=RiskHistogram.new,"out.data"=new.data))



  ####################################################
  ## Creat survival curves for patients of interest ##
  ####################################################

  if(!is.null(surv.print)){
    mut <- new.data[surv.print,]
    colnames(mut)[ncol(mut)] <- "RiskScore"

    allSurvs <- data.frame(nrow= 5)

    for(j in 1:nrow(mut)){
      survival.probs <- as.data.frame(matrix(nrow=6,ncol=15))
      rownames(survival.probs) <- c("Patient","Surv","Lower","Upper","Time","OncoRiskScore")
      surv.temp <- survfit(riskRefit, newdata = mut[j,])
      for(i in 1:ncol(survival.probs)){
        survival.probs[,i] <- try(c(rownames(mut)[j],as.numeric(summary(surv.temp, times = (i*3-3))$surv),
                                    round(summary(surv.temp, times = (i*3-3))$lower,digits=2),
                                    round(summary(surv.temp, times = (i*3-3))$upper,digits=2),
                                    i*3-3,as.numeric(mut$RiskScore[j])),silent=T)
      }
      allSurvs <- cbind(allSurvs,survival.probs)
    }
    allSurvs <- allSurvs[,-1]

    a <- list(
      autotick = FALSE,
      dtick = 6,
      tickcolor = toRGB("black")
    )

    t.survival.probs <- as.data.frame(t(allSurvs))
    for(k in 2:ncol(t.survival.probs)){
      t.survival.probs[,k] <- as.numeric(as.character(t.survival.probs[,k]))
    }
    y <- list(
      title = "Survival Probability"
    )
    IndSurvKM <- plot_ly(t.survival.probs, x = ~Time, y = ~Surv, name = ~Patient, type = 'scatter',
                         mode = 'lines+markers',hoverinfo="hovertext",color = ~Patient,
                         hovertext = ~paste("Genetic Risk Score :",round(OncoRiskScore,digits=3))
    ) %>% layout(yaxis = y,xaxis = ~a) %>%
      layout(xaxis = list(title = paste0("Time (Months)"), showgrid = TRUE),showlegend = FALSE) %>%
      add_ribbons(data = t.survival.probs,
                  ymin = ~Lower,
                  ymax = ~Upper,
                  line = list(color = 'rgba(7, 164, 181, 0.05)'),
                  fillcolor = 'rgba(7, 164, 181, 0.2)',
                  name = "Confidence Interval") #%>% layout(showlegend = FALSE)

  }
  else{IndSurvKM = NULL}
  return(list("data.out" = new.data,"RiskHist"=RiskHistogram.new,"IncKM" = IndSurvKM,"survivalEst"=t.survival.probs))
}






