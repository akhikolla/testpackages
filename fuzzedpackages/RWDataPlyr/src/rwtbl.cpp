//============================================================================
// Name        : rwdataplyr.cpp
// Author      : Alan Butler
// Version     :
// Copyright   : CC0
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

std::vector< std::vector<std::string> > parse_rdf_meta(std::vector<std::string> rdf) {
	std::vector< std::vector<std::string> > meta;
	std::vector<std::string> row (2);
	std::string line, token;
	int i = 0;
	size_t pos;

	while (rdf.at(i) != "END_PACKAGE_PREAMBLE") {
		line = rdf.at(i);
		i = i + 1;
		pos = line.find(":");

		if (pos != std::string::npos) {
			token = line.substr(0, pos);
			row.at(0) = token;
			token = line.substr(pos + 1, line.size() - 1);
			row.at(1) = token;

			meta.push_back(row);
		}
	}

	return(meta);
}

std::vector<std::string> parse_line(std::string line) {
	size_t pos;
	std::string token;
	std::vector<std::string> row;

	pos = line.find(":");

	if (pos != std::string::npos) {
		token = line.substr(0, pos);
		row.push_back(token);

		token = line.substr(pos + 1, line.size() - 1);
		// remove leading whitespace
		pos = token.find_first_not_of(" \t");
		token = token.substr(pos, token.size());
		row.push_back(token);
	}

	return(row);
}

std::vector< std::vector<std::string> > build_table(std::vector<std::string> vals, std::vector<std::string> timesteps) {
	std::vector < std::vector<std::string> > tmp;
	tmp.push_back(timesteps);
	
	if (vals.size() == 1) {
	  // then its a scalar slot
	  std::vector<std::string> v2(timesteps.size(), vals.at(0));
	  tmp.push_back(v2);
	} else {
	  tmp.push_back(vals);
	}
	
	return(tmp);
}

std::vector<std::string> get_year_month(std::string timestep) {
	std::string year;
  int month;
	std::vector<std::string> ym(2, "-99");
	std::string months[] = {"January","February","March","April","May","June",
    "July","August","September","October","November","December"};
	
	size_t pos, pos2;

	pos = timestep.find("-");
	if (pos != std::string::npos) {
		pos2 = timestep.find("-", pos + 1);

		if (pos2 != std::string::npos) {
			year = timestep.substr(0, 4);
			month = std::stoi(timestep.substr(pos + 1, pos2 - pos - 1));
			ym.at(0) = year;
			ym.at(1) = months[month - 1];
		}
	}

	return ym;
}

int get_n_runs(std::vector< std::vector<std::string> > meta) {
	// find the number_of_runs key word, get its value, and convert that to a size_t/int
	size_t i = 0;
	int val;
	while ((meta.at(i).at(0) != "number_of_runs") && (i < meta.size())) {
		i++;
	}

	val = std::stoi(meta.at(i).at(1));

	return val;
}

std::vector< std::vector<std::string> > parse_rdf(std::vector<std::string> rdf, int num_traces) {

	std::vector< std::vector<std::string> > atts, table, tmp_table;
	std::vector<std::string> row;
	std::string line, token;
	size_t i = 0;
	int trace;
	
	// skip meta section
	while (rdf.at(i) != "END_PACKAGE_PREAMBLE") {
		i += 1;
	}
	i += 1; // increment past END_PACKAGE_PREAMBLE

	for (trace = 1; trace <= num_traces; trace++) {
		std::vector<std::string> timesteps, ym, year, month;
		std::string slot_set, rule_set;
		size_t num_time_steps = 0;
		// read through END_RUN_PREAMBLE and parse value pairs based on ":"
		while(rdf.at(i) != "END_RUN_PREAMBLE") {
			line = rdf.at(i);
			i += 1;
			row = parse_line(line);

			// get the time_steps value and store as num_time_steps
			if (row.at(0) == "time_steps") {
				num_time_steps = std::stoi(row.at(1));
			} else if (row.at(0) == "slot_set") {
				slot_set = row.at(1);
			} else if (row.at(0) == "rule_set") {
				rule_set = row.at(1);
			} else {
				atts.push_back(row);
			}
		}
		i += 1; //increment past END_RUN_PREAMBLE

		// read the next time_steps rows and store as variable `timesteps`
		//    *** is this still correct for scalar slots
		for (size_t tt = 0; tt < num_time_steps; tt++) {
			timesteps.push_back(rdf.at(i));
			ym = get_year_month(rdf.at(i));
			year.push_back(ym.at(0));
			month.push_back(ym.at(1));
			i += 1;
		}

		while (rdf.at(i) != "END_RUN") {
			std::vector<std::string> vals;
			std::string obj, slot, obj_slot, object_type, units;

			// 2. read through END_SLOT_PREAMBLE and parse value pairs based on ":"
			//    - store object_type, object_name, slot_name, and create ObjectSlot
			while (rdf.at(i) != "END_SLOT_PREAMBLE") {
				line = rdf.at(i);
				i += 1;
				row = parse_line(line);
				if (row.at(0) == "object_name") {
					obj_slot = row.at(1);
					obj = row.at(1);
				} else if (row.at(0) == "slot_name") {
					obj_slot = obj_slot + "." + row.at(1);
					slot = row.at(1);
				} else if (row.at(0) == "object_type") {
					object_type = row.at(1);
				}  else {
					atts.push_back(row);
				}
			}

			i += 1;
			// read the next two lines and parse value pairs based on ":"
			row = parse_line(rdf.at(i));
			if (row.at(0) == "units") {
				units = row.at(1);
			}
			// *** need to add else that errors out with unexpected row

			i += 1;
			atts.push_back(parse_line(rdf.at(i))); // should be scale
			i += 1;

			// read the next num_time_steps rows, or read until you hit END_SLOT, skipping END_COLUMN
			while (rdf.at(i) != "END_SLOT") {
				if (rdf.at(i) != "END_COLUMN") {
					vals.push_back(rdf.at(i));
				}
				i += 1;
			}

			// create table here and append to existing table
			tmp_table = build_table(vals, timesteps);

			// ** add in check to make sure tmp_table.at(0) and tmp_table.at(1) are the same size
			if (table.size() > 0) {
				std::string trace_str = std::to_string(trace);
				for (size_t k = 0; k < tmp_table.at(0).size(); k++) {
					table.at(0).push_back(tmp_table.at(0).at(k));
					table.at(1).push_back(tmp_table.at(1).at(k));
					table.at(2).push_back(year.at(k));
					table.at(3).push_back(month.at(k));
					table.at(4).push_back(obj);
					table.at(5).push_back(slot);
					table.at(6).push_back(obj_slot);
					table.at(7).push_back(trace_str);
					table.at(8).push_back(slot_set);
					table.at(9).push_back(rule_set);
					table.at(10).push_back(object_type);
					table.at(11).push_back(units);
					// ** add in scale
				}
			} else {
				table.push_back(tmp_table.at(0));
				table.push_back(tmp_table.at(1));
				table.push_back(year);
				table.push_back(month);
				size_t nn = tmp_table.at(0).size();
				std::vector<std::string> col1(nn, obj), col2(nn, slot), col3(nn, obj_slot),
						col4(nn, std::to_string(trace)), col5(nn, slot_set), col6(nn, rule_set),
						col7(nn, object_type), col8(nn, units);

				table.push_back(col1);
				table.push_back(col2);
				table.push_back(col3);
				table.push_back(col4);
				table.push_back(col5);
				table.push_back(col6);
				table.push_back(col7);
				table.push_back(col8);
			}

			i += 1;
			// restart at 2.
		}
		i += 1;
	}
	return(table);
}

// [[Rcpp::export]]
List rdf_to_rwtbl_cpp(std::vector<std::string> rdf,
                      std::vector<std::string> keep_cols,
                      String const scenario = NA_STRING,
                      bool add_ym = true) {
  std::vector< std::vector<std::string> > meta, rwtbl;
  int num_runs, del_col;
  DataFrame val;
  List val_lst(3);
  
  
  meta = parse_rdf_meta(rdf);
  num_runs = get_n_runs(meta);
  rwtbl = parse_rdf(rdf, num_runs);
  
  size_t nn = rwtbl.at(0).size();
  // should maybe compare all the sizes...
  StringVector v0(nn), v3(nn), v4(nn), v5(nn), v6(nn), v8(nn), v9(nn), v10(nn),
    v11(nn), col_names;
  NumericVector v1(nn), v2(nn);
  IntegerVector v7(nn);
  std::vector<int> row_names(nn);
  
  v0 = rwtbl.at(0);
  v3 = rwtbl.at(3);
  v4 = rwtbl.at(4);
  v5 = rwtbl.at(5);
  v6 = rwtbl.at(6);
  v8 = rwtbl.at(8);
  v9 = rwtbl.at(9);
  v10 = rwtbl.at(10);
  v11 = rwtbl.at(11);
  
  for (size_t i = 0; i < nn; i++) {
    v1.at(i) = std::stod(rwtbl.at(1).at(i));
    v2.at(i) = std::stod(rwtbl.at(2).at(i));
    v7.at(i) = std::stoi(rwtbl.at(7).at(i));
    row_names.at(i) = i + 1;
  }
  
  
  // add in the scenario column if it's specified
  if (scenario != NA_STRING) {
    StringVector scen(nn, scenario);
    keep_cols.push_back("Scenario");
    val = DataFrame::create(
      _["Timestep"] = v0, 
      _["Year"] = v2,
      _["Month"] = v3,
      _["ObjectName"] = v4,
      _["SlotName"] = v5,
      _["ObjectType"] = v10,
      _["ObjectSlot"] = v6,
      _["Value"] = v1,
      _["Unit"] = v11,
      _["TraceNumber"] = v7,
      _["RulesetFileName"] = v9,
      _["InputDMIName"] = v8,
      _["Scenario"] = scen,
      _["stringsAsFactors"] = false );
  } else {
    val = DataFrame::create(
      _["Timestep"] = v0, 
      _["Year"] = v2,
      _["Month"] = v3,
      _["ObjectName"] = v4,
      _["SlotName"] = v5,
      _["ObjectType"] = v10,
      _["ObjectSlot"] = v6,
      _["Value"] = v1,
      _["Unit"] = v11,
      _["TraceNumber"] = v7,
      _["RulesetFileName"] = v9,
      _["InputDMIName"] = v8,
      _["stringsAsFactors"] = false );
  }
  
  // remove columns the user does not want
  if (add_ym) {
    keep_cols.push_back("Year");
    keep_cols.push_back("Month");
  }
  col_names = val.names();
  std::vector<std::string>::iterator it;
  for (int i = 0; i < col_names.size(); i++) {
    std::string col1 = Rcpp::as<std::string>(col_names.at(i));
    it = std::find(keep_cols.begin(), keep_cols.end(), col1);
    
    if (it == keep_cols.end()) {
      // column name was not found in keep_cols, so we want to remove it
      String tmp = col_names.at(i);
      del_col = val.findName(tmp);
      val.erase(del_col);
    }
  }
  
  val.attr("class") = "data.frame";
  val.attr("row.names") = row_names;
  val.attr("mrm_config_name") = meta.at(0).at(1);
  val.attr("owner") = meta.at(1).at(1);
  val.attr("description") = meta.at(2).at(1);
  val.attr("create_date") = meta.at(3).at(1);
  val.attr("n_traces") = std::stoi(meta.at(4).at(1));
  
  return val;
}
