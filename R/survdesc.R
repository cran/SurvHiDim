#' @title High dimensional univariate cox proportional hazard analysis.
#'
#' @param m Starting column number form where study variables of high dimensional data will get selected.
#' @param n Ending column number till where study variables of high dimensional data will get selected.
#' @param survdur Column name of survival duration event, a string value. i.e. "os"
#' @param event Column name of survival event, a string value. i.e "death"
#' @param aic By default aic = FALSE, if aic = TRUE the function returns
#' @param data High dimensional data containing the survival, progression and genomic observations.
#' @description Given the dimension of variables and survival information the function performs uni-variate Cox PH.
#' @import survival
#' @export
#' @return A data set containing estimates for variables present in column m to n.
#' @examples
#' ##
#' data(hnscc)
#'survdesc(m=10,n=50,survdur="os",event="death",aic=TRUE,data=hnscc)
#' ##
survdesc <- function(m,n,survdur,event,aic=TRUE,data){

  OS<-survdur
  Event<-event
  data1 <- data
  data2 <- subset(data, select=c(get(OS),get(Event),m:n))
  if(sum(is.na(data2))==0){
    data<-data2
  }

  if(sum(is.na(data2))>=1){
    message("Error in data : variable column contains missing values")
  }

  m=3
  n=ncol(data)

  #if((nrow(data)/n-m+1)>=10){
  da<-matrix(nrow = 0, ncol = 3) #creating dummy matrix
  dap<-matrix(nrow = 0, ncol = 1)
  aicval <- matrix(nrow = 0, ncol = 1)

    for(i in m:n){
      coxmod <- coxph(Surv(get(OS), get(Event)) ~ data[,i], data = data)
      dsum <- summary(coxmod)
      cia <- AIC(coxmod)
      da <- rbind(da,round(dsum$conf.int[,-2],4))
      dap <- rbind(dap,round(dsum$coefficients[,5],4))
      aicval <- rbind(aicval,cia)
    }

  dat <- data.frame(da,dap,aicval)
  nm <- colnames(data[m:n])
  dan <- cbind(nm,dat)
  colnames(dan)<- c("Variables","HR","LCL","UCL","Pvalue","AIC")
  rownames(dan)<-NULL

  if(aic == TRUE){
  dan <-dan[order(dan$AIC),]
  }
  output <- dan
  return(output)
}
utils::globalVariables(c("Status","pvalue","AIC"))
