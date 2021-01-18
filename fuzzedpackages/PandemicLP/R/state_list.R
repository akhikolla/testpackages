#' List of states available in the Covid-19 database
#'
#' @description This function provides the state abbreviations to be used in the \code{\link{load_covid}} function
#' and the corresponding state names. Only brazilian states are currently available.
#'
#' @seealso \code{\link{load_covid}} and \code{\link{country_list}}.
#'
#' @export
#'

state_list <- function(){

  state_abb <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA","MG","MS","MT",
                 "PA","PB","PE","PI","PR","RJ","RN","RO","RR","RS","SC","SE","SP","TO")
  state <- c("Acre","Alagoas","Amazonas","Amapa","Bahia","Ceara","Distrito Federal", "Espirito Santo","Goias",
             "Maranhao","Minas Gerais","Mato Grosso do Sul","Mato Grosso","Para","Paraiba","Pernambuco","Piaui",
             "Parana","Rio de Janeiro","Rio Grande do Norte","Rondonia","Roraima","Rio Grande do Sul",
             "Santa Catarina", "Sergipe", "Sao Paulo","Tocantins")

  state_list <- data.frame(state_abb, state)
  state_list

}
