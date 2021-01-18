#' Add date/time column to data.frame
#' 
#' Takes the `Year`, `Month`, `Day` and, if available, `Hour`, `Minute` and `Second` columns of a data.frame and uses them to produce a `Date` column that uses R's standard data/time format.
#' 
#' Converts YEARMODA to R date
#' 
#' @param df Data frame containing columns `Year`, `Month`, `Day` and - optionally - `Hour`, `Minute` and/or `Second`
#' @return data.frame consisting of the df input and a new column `Date`
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' 
#' add_date(KA_weather)
#' add_date(Winters_hours_gaps)
#' 
#'  
#' @export add_date
add_date<-function(df)
{
  if(!("Year" %in% colnames(df) & "Month" %in% colnames(df) & "Day" %in% colnames(df))) stop("Required input column 'Year', 'Month' and/or 'Day' is missing.")
  if(!"Hour" %in% colnames(df)) 
      df[,"Date"]<-ISOdate(df$Year,df$Month,df$Day) else
        if(!"Minute" %in% colnames(df)) 
          df[,"Date"]<-ISOdate(df$Year,df$Month,df$Day,df$Hour)  else  
            if(!"Second" %in% colnames(df)) 
              df[,"Date"]<-ISOdate(df$Year,df$Month,df$Day,df$Hour,df$Minute)  else  
                df[,"Date"]<-ISOdate(df$Year,df$Month,df$Day,df$Hour,df$Minute,df$Second) 
  
return(df) 
}
