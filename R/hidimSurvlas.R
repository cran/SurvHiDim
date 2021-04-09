#' @title hidimSurvlas
#'
#' Two step filteration without bonferroni correction
#'
#' @description  Survival analysis on high dimensional data using LASSO technique.
#'
#' @details hidimSurvlas function first creates a new column 'Status' in the input
#' data set and assigns values 0, 1, 2 to each rows. It assigns 0 (for progression = 1 & death(event) = 0) is
#' or (when progression = 0 & death(event) = 0. It assigns 1 (for progression = 1 &
#' death = 1, whereas assigns 2 (for progression = 0 and death = 1).
#'
#' Further, it creates two data sets, one data set named 'deathdata' which includes
#' subjects with status 0 and 1 and applies Cox PH on it. Another data is named as
#' 'compdata' which includes subjects with status 0 and 2, then applies Cox PH after
#' substituting 2 by 1. Then it filters out study variables having P-value < siglevel(significance level taken as input from user) from both
#' subset data. Secondly, it merges the commom sigificant variables from both data and creates a
#' new data frame which contains columns, 'ID','OS','Death','PFS','Prog','Status' and observations
#' of common significant variables (which are supposed to be leading to death given they leads to progression of cancer
#' as well as accounts for competing risks)
#' Further, it lists the common variables names and correspond in results in .csv format by default in user's current
#' working directory.
#'
#' On obtained commondata 'hidimSurvlas' fits LASSO and reduces the dimensions of study variables to handful significant variables.
#'
#' hidimSurvlas(m,n,siglevel,data)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) OS column must be named as 'OS'.
#'
#' 3) Death status/event column should be named as 'Death'.
#'
#' 4) Progression Free Survival column should be named as 'PFS'.
#'
#' 5) Progression event column should be named as 'Prog'.
#'
#' deathdata - A data frame with status column includes only those rows/subjects for which death/event was observed or not, given
#' progression was observed or not.
#'
#' compdata - A data frame with status column which includes those rows/subjects who died given progression was observed
#'
#' data1variables - List of variables/genes from deathdata
#'
#' data2variables - List of variables/genes from compdata
#'
#' siginificantpvalueA - A data frame with estimate values, HR, Pvalue, etc. of significant variables from deathdata.
#' siginificantpvalueB - A data frame with estimate values, HR, Pvalue, etc. of significant variables from compdata.
#'
#' commongenes - A data frame consisting observations on common significant study variables.
#'
#' cvar - List of common significant study variables.
#'
#' commondata - A final data out consisting survival information ans observations on common significant study variables.
#'
#' By default the function stores the output in .csv forms in current directory of user.
#'
#' @param m Starting column number form where study variables of high dimensional data will get selected.
#' @param n Ending column number till where study variables of high dimensional data will get selected.
#' @param data High dimensional data containing the survival, progression and genomic observations.
#' @param siglevel Level of significance pre-determined by the user
#' @param ID Column name of subject ID, a string value. i.e. "id"
#' @param OS Column name of survival duration event, a string value. i.e. "os"
#' @param Death Column name of survival event, a string value. i.e "death"
#' @param PFS Column name of progression free survival duration, a string value. i.e "pfs"
#' @param Prog Column name of progression event, a string value. i.e "prog"
#' @import survival
#' @import utils
#' @import glmnet
#' @import tidyverse
#' @importFrom  Rdpack reprompt
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#' @references Congdon, P. (2014). Applied bayesian modelling (Vol. 595). John Wiley & Sons.
#' @references Banerjee, S., Vishwakarma, G. K., & Bhattacharjee, A. (2019). Classification Algorithm for
#' High Dimensional Protein Markers in Time-course Data. arXiv preprint arXiv:1907.12853.
#' @return List of significant genes
#' @examples
#' ##
#' data(hnscc)
#' hidimSurvlas(7,105,0.05,ID="id",OS="os",Death="death",PFS="pfs",Prog="prog",hnscc)
#' ##
#' @export
#' @author Atanu Bhattacharjee and Akash Pawar
hidimSurvlas <- function(m,n,siglevel,ID,OS,Death,PFS,Prog,data){



  data$Status[data[,Prog] == 0 & data[,Death] == 1] <- 2
  data$Status[data[,Prog] == 0 & data[,Death] == 0] <- 0
  data$Status[data[,Prog] == 1 & data[,Death] == 1] <- 1
  data$Status[data[,Prog] == 1 & data[,Death] == 0] <- 0

  data1 <- data
  data2 <- data

  data2$Status[data2$Status == 1] <- 3 #death1 = 1 in mydata2 is assigned "3"
  data2$Status[data2$Status == 2] <- 1 #death1 = 2 in mydata2 is assigned "1"

  deathdata <- subset(data1, Status != 2) #deleting "2" row from mydata1
  compdata <- subset(data2, Status != 3) #deleting "3" row from mydata2

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
  nm <- colnames(data[m:n])
  dan <- cbind(nm,dat)
  dcn <- cbind(nm,dct)
  colnames(dan)<- c("Variables","HR","LCL","UCL","pvalue")
  colnames(dcn)<- c("Variables","HR","LCL","UCL","pvalue")
  dan <- subset(dan, pvalue < siglevel)
  dcn <- subset(dcn, pvalue < siglevel)


  significantpvalueA <- dan
  significantpvalueB <- dcn

  commongenes <- merge(significantpvalueA, significantpvalueB, by ="Variables")

  data1variables <- significantpvalueA$Variables
  data2variables <- significantpvalueB$Variables


  cvar <- commongenes$Variables
  Status ="Status"
  commondata1 <- data.frame(Idname<-data[ID],OverallSurvival<-data[OS],DeathEvent<-data[Death],ProgFreeSurv<-data[PFS],
                            ProgEvent<-data[Prog],StatusEvent<-data1[Status])
  commondata2 <- subset(data, select = c(cvar))

  commondata <-cbind(commondata1,commondata2)


   if(length(cvar)>0){
    OS<-commondata[,OS]
    Death<-commondata[,Death]
    commondata2=subset(commondata,select=c(cvar))
    x <- commondata2 %>% data.matrix()
    cv.commondata2 =cv.glmnet(x,Surv(OS,Death),alpha = 1,family = "cox",maxit =1000)
    plot(cv.commondata2)
    log(cv.commondata2$lambda.min)
    LambdaMin = cv.commondata2$lambda.min
    LambdaMax = cv.commondata2$lambda.max
    est.coef = coef(cv.commondata2, s = cv.commondata2$lambda.min)
    p <- as.matrix(est.coef)
    colnames(p) <- c("LambdaMin Value")
    s <- p[p[,1]!= 0,]
    s <- as.matrix(s)
    colnames(s) <- c("LambdaMin Value")
    l <- data.frame(rownames(s),s[,1])
    signif_genes <- list(rownames(s))

    print(paste(" Number of Significant DEGs are - ", nrow(s)))
    #===============================================================================

    rep <- list()
    rep[1] <- "Sigificant Genes after applying LASSO are - "

    rep[2] <- signif_genes
    return(rep)
    #write.csv(s, "Lasso_output.csv",row.names = T)
  }else{print("Cannot perform LASSO as second filteration as no gnes were significnt using first filteration which Cox Proportional Hazards, use HiDimSurv() to get minimum significant genes or co-variates")
    }
}
utils::globalVariables(c("Status","pvalue"))

