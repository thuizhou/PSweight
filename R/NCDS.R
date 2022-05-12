#' Illustrative dataset for PSweight
#'
#' This is a real observational study with binary and multiple treatment groups to illustrate the utility of PSweight functions.
#'
#' @docType data
#'
#' @usage data(NCDS)
#'
#' @format A data frame with 3642 rows and 16 columns.
#'
#' @details An dataset from the National Child Development Survey (NCDS) of the United Kingdom (UK). This dataset is obtained through the CC0 waiver from the work by Battistin and Sianesi (2011).
#'  For illustration, missing entries are imputed once with Multiple Imputation by Chained Equations (MICE).
#'
#' This data frame contains the following columns:
#' \itemize{
#' \item white: self-identified as white.
#'
#' \item wage: gross hourly wage in pound in log scale.
#'
#' \item Dany: whether received any education before.
#'
#' \item Dmult: three levels of educational attainment.
#'
#' \item maemp: employment status of mother.
#'
#' \item scht: school type.
#'
#' \item qmab: math score at age 7.
#'
#' \item qmab2: math score at age 11.
#'
#' \item qvab: reading score at age 7.
#'
#' \item qvab2: reading score at age 11.
#'
#' \item paed_u: father's years of education.
#'
#' \item maed_u: mother's years of education.
#'
#' \item age_pa: age of father.
#'
#' \item age_ma: age of mother.
#'
#' \item sub_u: number of siblings.
#'
#' \item wagebin: dichotomized wage obtained with a cutoff of average hourly wage.
#' }
#'
#' @references
#' https://cls.ucl.ac.uk/cls-studies/1958-national-child-development-study/
#'
#' Battistin E, Sianesi B. (2011). Misclassified Treatment Status and Treatment Effects: an Application
#' to Returns to Education in the United Kindom. Review of Economics and Statistics, 93(2), 495-509.
#'
#' Battistin E, Sianesi B. (2012). Replication data for: Misclassified Treatment Status and Treatment Effects: an Application
#' to Returns to Education in the United Kindom. URL https://doi.org.10.7910/DVN/EPCYUL.
#'
#' @keywords dataset
#'
#' @examples
#' data("NCDS")

"NCDS"
