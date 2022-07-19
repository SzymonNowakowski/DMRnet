#' miete dataset
#'
#' The miete data contains the rent index for Munich in 2003.
#'
#' @docType data
#'
#' @usage data(miete)
#'
#' @format A data frame with 2053 observations on the following 12 variables.
#' \describe{
#'   \item{rent}{Rent in euros.}
#'   \item{bathextra}{Special furniture in bathroom, yes = 1, no = 0.}
#'   \item{tiles}{Bathroom with tiles, yes = 0, no = 1.}
#'   \item{area}{Municipality.}
#'   \item{kitchen}{Upmarket kitchen, yes = 1, no = 0.}
#'   \item{rooms}{Number of rooms.}
#'   \item{best}{Best address, yes = 1, no = 0.}
#'   \item{good}{Good address, yes = 1, no =0.}
#'   \item{warm}{Warm water, yes = 0, no = 1.}
#'   \item{central}{Central heating, yes = 0, no = 1.}
#'   \item{year}{Year of construction.}
#'   \item{size}{Living space in square meter.}
#' }
#' @keywords datasets
#'
#' @references Fahrmeir, L., Künstler, R., Pigeot, I., Tutz, G. (2004) Statistik: der Weg zur Datenanalyse. 5. Auflage, Berlin: Springer-Verlag.
#'
#' @examples
#' data(miete)
#' summary(miete)
"miete"
