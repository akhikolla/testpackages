sail = function( df , code , fullData = TRUE , rowname = "_rowname_" , stringsAsFactors = FALSE ){
	colnames_df = colnames(df)
	rowname_added_temporarily = F

	if( rowname == F ){
	}else{
		rowname_within_ori_df = ( rowname %in% colnames_df )
		if( ! rowname_within_ori_df ){
			if( !is.null(nrow(df)) && ( nrow(df) == length(row.names(df)))){
				df[,rowname] = row.names(df)
				rowname_added_temporarily = T
			}else{
				cat("NOTE: `nrow(df)` returns inappropriate size: " , ifelse( is.null(nrow(df)), "NULL" , nrow(df) ) , "\n" , sep="" )
				cat("NOTE: rowname parameter is ignored. \n")
				# rowname needs not be added as a new column
			}
		}else{
			# rowname needs not be added as a new column
		}
	}
	# Add row number column for internal use
	colname_n = "_n_"
	if( colname_n %in% colnames(df) ){
		cat("NOTE: column '_n_' is replaced with row number for internal use, and will be deleted when execution is finished.\n")
	}
	df[,colname_n] = seq(1, nrow(df))

	# Add _discard_ column for internal use
	colname_discard = "_discard_"
	if( colname_discard %in% colnames(df) ){
		cat("NOTE: column '_discard_' is replaced with 0/1 for internal use, and will be deleted when execution is finished.\n")
	}
	df[,colname_discard] = as.integer( rep(0, nrow(df)))

	# If some columns have the same name, left most columns are used
	ori_colnames = colnames(df)
	ori_unique_colnames = unique(ori_colnames)

	positions_used_for_each_colname = sapply(ori_unique_colnames, function(unique_name){
		matched_positions = (ori_colnames %in% unique_name)
		if(sum(matched_positions) >= 2){
			cat(sprintf("Note: The original dataset has duplicated variable name( %s ). The leftmost column is used/updated.\n", unique_name))
		}
		min(which(matched_positions==TRUE ))
	}, USE.NAMES = TRUE )

	df_wo_duplicated_colnames = df[ positions_used_for_each_colname ]

	result = .data_sailr_cpp_execute( code, df_wo_duplicated_colnames)

	if( ("DataSailr_NewOrder" %in% names(attributes(result))) && attr(result, "DataSailr_NewOrder") == TRUE ){
		if(fullData == TRUE){
			order_vec = attr(result, "DataSailr_NewOrderVector")
			df = df[order_vec,]
		}
		attr(result, "DataSailr_NewOrder") = NULL
		attr(result, "DataSailr_NewOrderVector") = NULL
	}

	if( ncol(result) == 0 ){
		# If no assignments were executed in script. 
		assign_occured_in_script = FALSE
	}else{
		assign_occured_in_script = TRUE
	}

	if(stringsAsFactors == TRUE ){
		result_df = data.frame( lapply(result, function(x) if (is.factor(x)) as.character(x) else {x} ), stringsAsFactors = TRUE )  # Deal strings as factors
	}else{
		result_df = result # Deal strings as chracter vectors (Default action)
	}

	if( rowname_added_temporarily  ){
		df = df[ !( colnames(df) %in%  rowname )]
	}

	cols_for_update = colnames(result_df) %in% ori_unique_colnames #logical
	colnames_for_update = colnames(result_df)[cols_for_update] #character
	cols_for_addition =  !cols_for_update #logical
	if(fullData == T){
		if(assign_occured_in_script){
			# update original columns
			lapply(colnames_for_update, function(colname_for_update){
				pos_to_update = positions_used_for_each_colname[colname_for_update]
				df[pos_to_update] <<- result_df[colname_for_update]
			})
			# add new columns
			result_df = cbind(df[, -which(names(df) %in% c("_n_", "_discard_"))] , result_df[cols_for_addition])
		}else{
			# If no assignments were executed in script, use the original dataframe with updated row orders.
			# Also, _n_ and _discard_ are removed.
			result_df = df[, -which(names(df) %in% c("_n_", "_discard_"))]
		}
	}

	return(result_df)
}


author = function(){
  print("DataSailr is actively developed by Toshi Umehara (@niceume).")
}
