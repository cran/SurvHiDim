#' @title hnscc
#'
#' High dimensional genomic data on head and neck cancer
#'
#' @description High dimensional breast cancer gene expression data
#' @usage hnscc
#' @format A dataframe with 565 rows and 104 variables
#' \describe{
#' \item{id}{ID of subjects}
#' \item{leftcensor}{Initial censoring time}
#' \item{death}{Survival event}
#' \item{os}{Duration of overall survival}
#' \item{pfs}{Duration of progression free survival}
#' \item{prog}{Progression event}
#' \item{GJB1,...,HMGCS2}{High dimensional covariates}}
#' @examples
#' \dontrun{
#' data(hnscc)
#' }
"hnscc"

