# checks the input object data

stratEst.check.data <- function( data ){

  if( "data.frame" %in% class(data) == F ){
    stop("stratEst error: The input object 'data' must be an object of class 'stratEst.data'.")
  }

  # check mandatory data frame variables
  if( is.null(data$id) ) {
    stop("stratEst error: Data does not contain the variable 'id'.")
  }
  if( is.null(data$game) ) {
    stop("stratEst error: Data does not contain the variable 'game'.")
  }
  if( is.null(data$period) ) {
    stop("stratEst error: Data does not contain the variable 'period'.")
  }
  input.is.null = F
  if( is.null(data$input) ) {
    input.is.null = T
    data$input <- 1
  }
  if( is.null(data$choice) ) {
    stop("stratEst error: Data does not contain the variable 'choice'.")
  }
  else{

  }

  id <- data$id
  game <- data$game
  period <- data$period
  input <- data$input
  output <- data$choice

  # id
  id_is_factor = FALSE
  if( "factor" %in% class(id) ){
    id_factor <- id
    id_is_factor = TRUE
  }
  id <- as.numeric(id)
  if( any( is.na( id ) ) ){
    stop("stratEst error: The variable 'id' in data cannot contain NA values.");
  }

  # input
  if( "factor" %in% class(input) ){
    input_factor <- input
  }
  else{
    input_factor <- as.factor( input )
  }
  input <- match(input_factor,sort(unique(input_factor)))
  input[is.na(input)] <- 0
  input <- as.numeric(input)
  levels_input <- levels(input_factor)

  # output
  if( "factor" %in% class(output) ){
    output_factor <- output
  }
  else{
    output_factor <- as.factor( output )
  }
  output <- match(output_factor,sort(unique(output_factor)))
  levels_output <- levels( output_factor )
  if( any( is.na( output ) ) ){
    stop("stratEst error: The variable 'choice' in data cannot contain NA values.");
  }

  stratEst.check.data.return <- list( "data" = data , "id" = id , "game" = game , "period" = period , "input" = input , "output" = output, "input.factor" = input_factor , "output.factor" = output_factor , "levels.input" = levels_input , "levels.output" = levels_output, "input.is.null" = input.is.null )

  return(stratEst.check.data.return)
}
