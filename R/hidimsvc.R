#' @title hidimsvc
#' Survival analysis on high dimensional by creating batches of variables
#'
#' @description Survival analysis on high dimensional data by creating batches of covariates
#'
#' @details hidimsvc function fits Univarate Cox Proportinal Hazard models by considering each variables at a time. Then it filters out study variables having P-value < siglevel(significance level taken as input from user).
#' Once by survival and survival eevent and another by progression and progression events.
#' Secondly, it merges the commom sigificant variables from both OS and PFS analysis and creates a
#' new data frame which contains columns, 'ID','OS','Death','PFS','Prog','Status' and observations
#' of common significant variables (which are supposed to be leading to death given they leads to progression of cancer
#' as well as accounts for competing risks)
#' Further, it lists the common variables names and outs corresponding results in .csv format by default in user's current
#' working directory.
#'
#' It works similary to HiDimSurv unlike it creates batches of decided study vriables by user to make the analysis less time consuming.
#'
#' hidimsvc(m1,m2,batchsize,siglevel,data),
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) OS column must be named as 'OS'.
#'
#' 3) Death status/event column should be named as 'Death'.
#'
#' 4) Progression Fress Survival column should be named as 'PFS'.
#'
#' 5) Progression event column should be named as 'Prog'.
#'
#' OSDeathcoeff - A data frame containing HR estimates and p-values for study variables on fitting univariate CoxPh on OS and Survival event.
#'
#' PFSProgcoeff - A data frame containg HR estimates and p-values for study variables on fitting univariate CoxPh on PFS and Progression event.
#'
#' namevect - List of all the study variable names.
#'
#' significantOSDeathgenes - A data frame containing HR estimates and p-values for significant study variables.
#'
#' significantPFSProggenes - A data frame containing HR estimates and p-values for significant study variables
#'
#' commongenes - A data frame containing estimated values of significant study variables found common from significant study variables
#' on fitting CoxPh on survival and progression times and events.
#'
#' odnames - List of significant variables on fitting CoxPh using survival and survival event.
#'
#' ppnames - List of significant variables on fitting CoxPh using progression and progression event.
#'
#' cvar - List of common significant study variables on fitting CoxPh on survival and progression, times and events.
#'
#' commondata - A data out which contains the clinical observations and observations on commongenes variables.
#'
#' @param m Starting column number form where study variables of high dimensional data will get selected.
#' @param n Ending column number till where study variables of high dimensional data will get selected.
#' @param batchsize Number of variables to be consider at time while running function (maximum batch size should not be greater than one third of the total number of high dimensional variables)
#' @param siglevel Level of significance pre-determined by the user
#' @param ID Column name of subject ID, a string value. i.e. "id"
#' @param OS Column name of survival duration event, a string value. i.e. "os"
#' @param Death Column name of survival event, a string value. i.e "death"
#' @param PFS Column name of progression free survival duration, a string value. i.e "pfs"
#' @param Prog Column name of progression event, a string value. i.e "prog"
#' @param data High dimensional data containing the survival, progression and genomic observations.
#' @import survival
#' @import utils
#' @importFrom  Rdpack reprompt
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#' @references Congdon, P. (2014). Applied bayesian modelling (Vol. 595). John Wiley & Sons.
#' @references Banerjee, S., Vishwakarma, G. K., & Bhattacharjee, A. (2019). Classification Algorithm for
#' High Dimensional Protein Markers in Time-course Data. arXiv preprint arXiv:1907.12853.
#' @return Estimate values of significant variables/DEGs on considering Death with Progression
#' @return Estimate values of significant variables/DEGs on considering Death without Progression
#' @return List of variable/DEGs considering Death with Progression
#' @return List of variable/DEGs considering Death without Progression
#' @return Estimates data for the DEGs/Variables found common between significant DEGs from data having death due to progression and data showing death without progression"
#' @return List of variable/DEGs found common between significant DEGs from data having death due to progression and data showing death without progression"
#' @examples
#' ##
#' data(hnscc)
#' hidimsvc(7,105,5,0.05,ID="id",OS="os",Death="death",PFS="pfs",Prog="prog",hnscc)
#' ##
#' @export
#' @author Atanu Bhattacharjee and Akash Pawar

hidimsvc <- function(m,n,batchsize,siglevel,ID,OS,Death,PFS,Prog,data){

  m=m;
  m2=n
  d <- data
  data1 <- d[,c(1:m2)]

  k = batchsize
  n=m+(k-1)

  nbatch <- (ncol(data1) - (m-1))/k
  da<-matrix(nrow = 0, ncol = 3) #creating dummy matrix
  dc<-matrix(nrow = 0, ncol = 3)
  dap<-matrix(nrow = 0, ncol = 1)
  dcp<-matrix(nrow = 0, ncol = 1)
  for(i in m:n){

    b = data[,i]
    c = data[,i]
    coxdeath <- coxph(Surv(get(OS), get(Death)) ~ b, data = data)
    coxcomp <- coxph(Surv(get(PFS), get(Prog)) ~ c, data = data)
    dsum <- summary(coxdeath)
    csum <- summary(coxcomp)

    da <- rbind(da,dsum$conf.int[,-2])
    dc <- rbind(dc,csum$conf.int[,-2])

    dap <- rbind(dap,round(dsum$coefficients[,5],4))
    dcp <- rbind(dcp,round(csum$coefficients[,5],4))
  }
  dat <- data.frame(da,dap)
  dct <- data.frame(dc,dcp)
  rownames(dat)=NULL;rownames(dct)=NULL
  nm <- colnames(data1[m:n])
  dan <- cbind(nm,dat)
  dcn <- cbind(nm,dct)
  colnames(dan)<- c("Variables","HR","LCL","UCL","pvalue")
  colnames(dcn)<- c("Variables","HR","LCL","UCL","pvalue")
  dan <- subset(dan, pvalue < siglevel)
  dcn <- subset(dcn, pvalue < siglevel)


  l1 = list()
  l2 = list()

  int <- seq(1,nbatch,1)
  int1 <- seq(1,nbatch,1)
  init <- length(int)
  p<-(m+(k*int[length(int)])):((m-1)+k*(int[length(int)]+1))

  if(p[length(p)]==ncol(data1)){
    int = int
  } else {int = int[1:int[length(int)-1]]}

  for(j in int){
    if(((ncol(data1)-(m-1))-(k*j)) >= k){
      d1<-matrix(nrow = 0, ncol = 3)
      d2<-matrix(nrow = 0, ncol = 3)
      cn1<-matrix(nrow = 0, ncol = 1)
      cn2<-matrix(nrow = 0, ncol = 1)
      d1p<-matrix(nrow = 0, ncol = 1)
      d2p<-matrix(nrow = 0, ncol = 1)
      for(i in (m+(k*j)):((m-1)+k*(j+1))){

        b = data[,i]
        c = data[,i]

        coxdeath <- coxph(Surv(get(OS), get(Death)) ~ b, data = data)
        coxcomp <- coxph(Surv(get(PFS), get(Prog)) ~ c, data = data)
        dsum <- summary(coxdeath)
        csum <- summary(coxcomp)
        d1 <- rbind(d1,dsum$conf.int[,-2])
        cn1 <- rbind(cn1,colnames(data1[i]))
        d2 <- rbind(d2,dsum$conf.int[,-2])
        cn2 <- rbind(cn2,colnames(data1[i]))

        d1p <- rbind(d1p,dsum$coefficients[,5])
        d2p <- rbind(d2p,csum$coefficients[,5])

      }
      dp1 <- data.frame(cn1,d1,d1p)
      dp2 <- data.frame(cn2,d2,d2p)
      rownames(dp1)=NULL;rownames(dp2)=NULL
      #dp3 <- cbind(cn1,dp1)
      dp3 <- dp1
      dp4 <- dp2
      #dp4 <- cbind(cn2,dp2)
      colnames(dp3)<- c("Variables","HR","LCL","UCL","pvalue")
      colnames(dp4)<- c("Variables","HR","LCL","UCL","pvalue")
      dp3 <- subset(dp3, pvalue < siglevel)
      dp4 <- subset(dp4, pvalue < siglevel)
      l1[[j]]<-dp3
      l2[[j]]<-dp4
    }
    else {
      d3<-matrix(nrow = 0, ncol = 3)
      d4<-matrix(nrow = 0, ncol = 3)
      d3pv<-matrix(nrow = 0,ncol =1)
      d4pv<-matrix(nrow = 0,ncol =1)

      u = (m+(k*j))
      v =((m+(k*j))+(ncol(data1)-(m-1))-(k*j)-1)
      for(i in u:v){
        b = data[,i]
        c = data[,i]
        coxdeath <- coxph(Surv(get(OS), get(Death)) ~ b, data = data)
        coxcomp <- coxph(Surv(get(PFS), get(Prog)) ~ c, data = data)
        dsum <- summary(coxdeath)
        csum <- summary(coxcomp)

        d3 <- rbind(d3,dsum$conf.int[,-2])
        d4 <- rbind(d4,dsum$conf.int[,-2])
        d3pv <- rbind(d3pv,dsum$coefficients[,5])
        d4pv <- rbind(d4pv,dsum$coefficients[,5])
      }
      nm3 <- colnames(data1[(m+(k*j)):((m+(k*j))+(ncol(data1)-(m-1))-(k*j)-1)])

      d3t <- data.frame(nm3,d3pv,d3)
      d4t <- data.frame(nm4,d4pv,d4)
      #d3n <- cbind(nm3,d3t)

      colnames(d3t)<- c("Variables","HR","LCL","UCL","pvalue")
      #d4n <- cbind(nm3,d4t)
      colnames(d4t)<- c("Variables","HR","LCL","UCL","pvalue")
      rownames(d3n)=NULL;rownames(d4n)=NULL

    }
  }

  fda <- matrix(nrow = 0, ncol = 5)
  fdc <- matrix(nrow = 0, ncol = 5)


  for(j in 1:(length(int)-1)){
    fda <- rbind(fda,l1[[j]])
    fdc <- rbind(fdc,l2[[j]])
  }
  d1t <- data.frame(fda)
  d2t <- data.frame(fdc)

  if(init != length(int)){
    OSDeathcoeff <- rbind(dan,d1t)
    PFSProgcoeff <- rbind(dcn,d2t)
  }else{
    OSDeathcoeff <- rbind(dan,d1t,d3n)
    PFSProgcoeff <- rbind(dcn,d2t,d4n)
  }

  GeneOSDeath <- function(){
    return(OSDeathcoeff)
  }
  GenePFSProg <- function(){
    return(PFSProgcoeff)
  }
  namevect <- names(data1)

  sigOSDeathgenes <- function(){
    return(significantOSDeathgenes)
  }
  sigPFSProggenes <- function(){
    return(significantPFSProggenes)
  }

  commongenes <- merge(OSDeathcoeff, PFSProgcoeff, by ="Variables")
  commonDeathProggenes <- function(){
    return(commongenes)
  }
  odnames <- OSDeathcoeff$Variables
  ppnames <- PFSProgcoeff$Variables

  OSDeathGenes <-function(){
    return(odnames)
  }
  PFSProgGenes <-function(){
    return(ppnames)
  }

  cvar <- commongenes$Variables
  commonOSPFSgenes <- function(){
    return(cvar)
  }

  commondata1 <- data.frame(ID<-data[ID],OS<-data[OS],Death<-data[Death],PFS<-data[PFS],
                            Prog<-data[Prog])
  commondata2 <- subset(data, select = c(cvar))

  commondata <-cbind(commondata1,commondata2)


  eventcompetdata <- function(){
    return(commondata)
  }
  print("Function executed properly")

  report1 <- list()
  report1[[1]] <- "Estimate values of significant variables/DEGs on considering Death
    with Progression"
  report1[[2]] <- OSDeathcoeff
  print(report1)

  report2 <- list()
  report2[[1]] <- "Estimate values of significant variables/DEGs considering Death without
    Progression data"
  report2[[2]] <- PFSProgcoeff
  print(report2)

  report3 <- list()
  report3[[1]] <- "List of variable/DEGs considering Death with Progression data"
  report3[[2]] <- odnames
  print(report3)

  report4 <- list()
  report4[[1]] <- "List of variable/DEGs considering Death without Progression data"
  report4[[2]] <- ppnames
  print(report4)

  #  report5 <- list()
  #  report5[[1]] <- "Estimates values for significant variables/DEGs from Deathdata/Death
  #    with Progression data"
  #  report5[[2]] <- significantOSDeathgenes
  #  print(report5)

  #  report6 <- list()
  #  report6[[1]] <- "Estimates for significant variables/DEGs from Compdata/Death without
  #    Progression data"
  #  report6[[2]] <- significantPFSProggenes
  #  print(report6)

  report5 <- list()
  report5[[1]] <- "Estimates data for the DEGs/Variables found common between significant
    DEGs from data having death due to progression and data showing death without progression"
  report5[[2]] <- commongenes
  print(report5)

  report6 <- list()
  report6[[1]] <- "List of variable/ DEGs found common between significant
    DEGs from data having death due to progression and data showing death without progression"
  report6[[2]] <- cvar
  print(report6)
  output<-list(report1,report2,report3,report4,report5,report6)
  return(output)
}
utils::globalVariables(c("pvalue","nm4","significantOSDeathgenes","significantPFSProggenes"))

