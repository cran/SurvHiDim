#' @title hidimSurvbon
#'
#' Uses bonferroni correction factor in survival analysis
#'
#' @description Applies HiDimSurv on high dimensional data using Bonferroni's correction criteria.
#'
#' @details hidimSurvbon function first creates a new column 'Status' in the input
#' data set and assigns values 0, 1, 2 to each rows. It assigns 0 (for progression = 1 & death(event) = 0) is
#' or (when progression = 0 & death(event) = 0. It assigns 1 (for progression = 1 &
#' death = 1, whereas assigns 2 (for progression = 0 and death = 1).
#'
#' Further, it creates two data sets, one data set named 'deathdata' which includes
#' subjects with status 0 and 1 and applies Cox PH on it. Another data is named as
#' 'compdata' which includes subjects with status 0 and 2, then applies Cox PH after
#' substituting 2 by 1. Then it filters out study variables having P-value < siglevel/no. of columns of high dimensional data(significance level taken as input from user) from both
#' subset data. Secondly, it merges the commom sigificant variables from both data and creates a
#' new data frame which contains columns, 'ID','OS','Death','PFS','Prog','Status' and observations
#' of common significant variables (which are supposed to be leading to death given they leads to progression of cancer
#' as well as accounts for competing risks)
#' Further, it lists the common variables names and correspond in results in .csv format by default in user's current
#' working directory.
#'
#'
#' hidimSurvbon(m,n,siglevel,data),
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
#' deathdata - A data frame with status column includes only those rows/subjects for which death/event was observed or not, given
#' progression was observed or not.
#'
#' compdata - A data frame with status column which includes those rows/subjects who died given progression was observed
#'
#' data1variables - list of variables/genes from deathdata
#'
#' data2variables - list of variables/genes from compdata
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
#' By default the fucntion stores the output in .csv forms in current directory of user.
#'
#' @param m Starting column number form where study variables of high dimensional data will get selected.
#' @param n Ending column number till where study variables of high dimensional data will get selected.
#' @param boncorr Level of significance on which bonferroni correction will be applied.
#' @param data High dimensional data containing the survival, progression and genomic observations.
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
#' @return List of variables found significant on OS and survival event
#' @return List of variables found significant on PFS and progression event
#' @return Estimates values for significant variables on OS and survival event
#' @return Estimates values for significant variables on PFS and progression event
#' @return Estimates data for the DEGs/Variables found common between significant DEGs from data having death due to progression and data showing death without progression
#' @return Estimates data for the DEGs/Variables found common between significant DEGs from data having death due to progression and data showing death without progression
#' @return List of variable/ DEGs found common between significantDEGs from data having death due to progression and data showing death without progression
#' @return Data with Survival outcomes and DEGs/Variable observations on each subject for DEGs found playing crucial role in death due to progression and without

#' @examples
#' ##
#' data(hnscc)
#' hidimSurvbon(7,105,0.05,ID="id",OS="os",Death="death",PFS="pfs",Prog="prog",hnscc)
#' ##
#' @export
#' @author Atanu Bhattacharjee and Akash Pawar
#' @seealso hidimSurvbonlas

hidimSurvbon <- function(m,n,boncorr,ID,OS,Death,PFS,Prog,data){


  siglevel <- boncorr


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
  dan <- subset(dan, pvalue < (siglevel/length(m:n)))
  dcn <- subset(dcn, pvalue < (siglevel/length(m:n)))


  significantpvalueA <- dan
  significantpvalueB <- dcn

  commongenes <- merge(significantpvalueA, significantpvalueB, by ="Variables")

  data1variables <- significantpvalueA$Variables
  data2variables <- significantpvalueB$Variables


  cvar <- commongenes$Variables
  Status ="Status"
  commondata1 <- data.frame(ID<-data[ID],OS<-data[OS],Death<-data[Death],PFS<-data[PFS],
                            Prog<-data[Prog],Status<-data1[Status])
  commondata2 <- subset(data, select = c(cvar))

  commondata <-cbind(commondata1,commondata2)


  report1 <- list()
  report1[[1]] <- "List of variables found significant on OS and survival event"
  report1[[2]] <- data1variables
  print(report1)

  report2 <- list()
  report2[[1]] <- "List of variables found significant on PFS and progression event"
  report2[[2]] <- data2variables
  print(report2)

  report3 <- list()
  report3[[1]] <- "Estimates values for significant variables on OS and survival event"
  report3[[2]] <- significantpvalueA
  print(report3)

  report4 <- list()
  report4[[1]] <- "Estimates values for significant variables on PFS and progression event"
  report4[[2]] <- significantpvalueB
  print(report4)

  report5 <- list()
  report5[[1]] <- "Estamates data for the DEGs/Variables found common between significant
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

  if(length(cvar)>=1){report7 <- list()
  report7[[1]] <- "Data with Survival outcomes and DEGs/Variable observations on each subject for
    DEGs found playing crucial role in death due to progression and without."
  report7[[2]] <- commondata
  print(report7)}


  if(length(cvar)<1){
    print("Using bonferroni criteria no variables were found to be significant, use HiDimSurv() function.")
  }

}
utils::globalVariables(c("Status","pvalue"))
