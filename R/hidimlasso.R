#' @title hidimLasso
#'
#' LASSO for high dimensional data
#'
#' @description Least Absolute Shrinkage and Selection Operator (LASSO) for High Dimensional Survival data.
#'
#' @details 'HiDimLasso' allows a user to apply LASSO function on the High Dimensional data and reduce the study
#' variables to handful number of co-variate which are observed impacting the survival outcomes.
#'
#' Column of Overall Survival must be named as 'OS' and the column defining the event
#' must be named as 'Event'.
#'
#' By default it stores the outcome data in user's current directory.
#' @param m Starting column number from where variables of high dimensional data will get selected.
#' @param n Ending column number till where variables of high dimensional data will get selected.
#' @param OS Column name of survival duration event, a string value. i.e. "os"
#' @param Death Column name of survival event, a string value. i.e "death"
#' @param data High dimensional data having survival duration, event and various covariates observations
#' @import survival
#' @import utils
#' @import glmnet
#' @importFrom stats coef
#' @author Atanu Bhattacharjee and Akash Pawar
#' @export
#' @return A list of variables selected by LASSO as predictor variables.
#' @examples
#' ##
#' data(hnscc)
#' hidimLasso(7,105,OS="os",Death="death",hnscc)
#' ##
#' @seealso hidimSurvlas hidimSurvbonlas
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#' @references Congdon, P. (2014). Applied bayesian modelling (Vol. 595). John Wiley & Sons.
#' @references Banerjee, S., Vishwakarma, G. K., & Bhattacharjee, A. (2019). Classification Algorithm for
#' High Dimensional Protein Markers in Time-course Data. arXiv preprint arXiv:1907.12853.
hidimLasso <- function(m,n,OS,Death,data){
  u=m;v=n;

  OS<-data[,OS]
  Death<-data[,Death]
  commondata2=subset(data,select=c(u:v))
  x <- as.matrix(commondata2)
  cv.commondata2 =cv.glmnet(x,Surv(OS,Death),alpha = 1,family = "cox",maxit =1000)

  log(cv.commondata2$lambda.min)
  LambdaMin = cv.commondata2$lambda.min
  LambdaMax = cv.commondata2$lambda.1se
  est.coef = coef(cv.commondata2, s = cv.commondata2$lambda.min)
  p <- as.matrix(est.coef)
  colnames(p) <- c("LambdaMin Value")
  s <- p[p[,1]!= 0,]
  s <- as.matrix(s)
  colnames(s) <- c("LambdaMin Value")
  l <- data.frame(rownames(s),s[,1])
  signif_genes <- list(rownames(s))

  if(nrow(l)==0){
    print("Chosen number of variables are not significant, increase the number of columns")
  } else{

  colnames(l)<-c("Sr","Lambda")
  Sr2<-seq(1,length(l$Lambda),by =1)
  plotdata<- data.frame(Sr2,l$Lambda)
  plot(plotdata)
  nrow(s)
  print(paste(" Number of Significant DEGs are - ", nrow(s)))
#===============================================================================
  plot(cv.commondata2)
  rep <- list()
  rep[1] <- "Sigificant Genes are - "
  rep[2] <- signif_genes
  return(rep)
  }
}

