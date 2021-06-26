#' @title High dimensional multivariate cox proportional hazard data analysis
#'
#' @param C1 Covar1
#' @param C2 Covar2
#' @param C3 Covar3
#' @param C4 Covar4
#' @param C5 Covar5
#' @param OS "Column/Variable name" consisting duration of survival.
#' @param event "Column/Variable name" consisting survival event.
#' @param data High dimensional data containing survival observations and covariates.
#' @description Given the dimensions of the variables and survival informations. The function
#' performs multivariate Cox PH by taking 5 variables at a time.
#'
#' @return Data set containing the survival estimates and Pvalue.
#' @import survival
#' @export
#'
#' @examples
#' ##
#' multicoxb(C1="GJB1",C2=NULL,C3="HPN",C4=NULL,C5=NULL,OS="os",event="death",data=hnscc)
multicoxb <- function(C1=NULL,C2=NULL,C3=NULL,C4=NULL,C5=NULL,OS,event,data){
  covar1<-C1
  covar2<-C2
  covar3<-C3
  covar4<-C4
  covar5<-C5
  Surv<- OS
  Event<-event
  data1 <- data
  multi <- c(covar1,covar2,covar3,covar4,covar5)

  data2a <- subset(data, select=c(get(Surv),get(Event)))
  data2b <- subset(data, select=c(multi))
  data2 <- data.frame(data2a,data2b)

  data <- data2

  if(length(multi)==5){
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]]+data[,multi[3]]+data[,multi[4]]+data[,multi[5]], data=data)
    sumr <- summary(model1)
    sumrcoeff<-round(sumr$coefficients[,c(2,5)],2)
    hr <- round(sumr$conf.int[,-2],2)

    p1 <- round(sumr$coefficients[,5][1],2)
    p2 <- round(sumr$coefficients[,5][2],2)
    p3 <- round(sumr$coefficients[,5][3],2)
    p4 <- round(sumr$coefficients[,5][4],2)
    p5 <- round(sumr$coefficients[,5][5],2)
    pv <- rbind(p1,p2,p3,p4,p5)

    aicval <- round(AIC(model1),4)
    cn <- multi
    fd <- cbind(cn,hr,pv,aicval)
    fadata <- data.frame(fd)
    rownames(fadata)<-NULL
    colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue","AIC")
  }
  if(length(multi)==4){
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]]+data[,multi[3]]+data[,multi[4]], data=data)
    sumr <- summary(model1)
    sumrcoeff<-round(sumr$coefficients[,c(2,5)],2)
    hr <- round(sumr$conf.int[,-2],2)

    p1 <- round(sumr$coefficients[,5][1],2)
    p2 <- round(sumr$coefficients[,5][2],2)
    p3 <- round(sumr$coefficients[,5][3],2)
    p4 <- round(sumr$coefficients[,5][4],2)
    pv <- rbind(p1,p2,p3,p4)
    aicval <- round(AIC(model1),4)
    cn <- multi
    fd <- cbind(cn,hr,pv,aicval)
    fadata <- data.frame(fd)
    rownames(fadata)<-NULL
    colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue","AIC")
  }
  if(length(multi)==3){
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]]+data[,multi[3]], data=data)
    sumr <- summary(model1)
    sumrcoeff<-round(sumr$coefficients[,c(2,5)],2)
    hr <- round(sumr$conf.int[,-2],2)

    p1 <- round(sumr$coefficients[,5][1],2)
    p2 <- round(sumr$coefficients[,5][2],2)
    p3 <- round(sumr$coefficients[,5][3],2)
    pv <- rbind(p1,p2,p3)
    aicval <- round(AIC(model1),4)
    cn <- multi
    fd <- cbind(cn,hr,pv,aicval)
    fadata <- data.frame(fd)
    rownames(fadata)<-NULL
    colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue","AIC")
  }
  if(length(multi)==2){
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]]+data[,multi[2]], data=data)
    sumr <- summary(model1)
    sumrcoeff<-round(sumr$coefficients[,c(2,5)],2)
    hr <- round(sumr$conf.int[,-2],2)
    aicval <- round(AIC(model1),4)
    p1 <- round(sumr$coefficients[,5][1],2)
    p2 <- round(sumr$coefficients[,5][2],2)
    pv <- rbind(p1,p2)
    cn <- multi
    fd <- cbind(cn,hr,pv,aicval)
    fadata <- data.frame(fd)
    rownames(fadata)<-NULL
    colnames(fadata)<-c("Variables","HR","LCL","UCL","Pvalue","AIC")
  }
  if(length(multi)==1){
    model1 <- coxph(Surv(get(Surv),get(Event)) ~ data[,multi[1]], data=data)
    sumr <- summary(model1)
    sumrcoeff<-round(sumr$coefficients[,c(2,5)],2)
    hr <- round(sumr$conf.int[,-2],2)
    aicval <- round(AIC(model1),4)
    p1 <- round(sumr$coefficients[,5][1],2)


    cn <- multi

    fadata <- data.frame(Variables = cn, HR = hr[1], LCL = hr[2], UCL = hr[3], Pvalue = p1, AIC = aicval)
    rownames(fadata)<-NULL
  }
  return(fadata)
}

utils::globalVariables(c("C1","C2","C3","C4","C5"))
