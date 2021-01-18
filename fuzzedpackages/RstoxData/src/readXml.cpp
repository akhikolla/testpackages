#include <iostream>
#include <fstream>
#include <vector>
#include <cstddef>
#include <cassert>
#include <map>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

#define PUGIXML_HEADER_ONLY
#define PUGIXML_NO_EXCEPTIONS

#include "pugixml/pugixml.hpp"

#define STRICT_R_HEADERS
#include "Rcpp.h"

void processNode(pugi::xml_node& node, const std::vector<const char*>& parentPrefix, std::map<std::string, std::vector<std::string> >& tableHeaders, std::map<std::string, int >& prefixLens, std::map<std::string, int>& levelCtrs, Rcpp::List& ret) {

	const char* root = node.name();

	// If we want to skip this node, go to sibling
	while (tableHeaders.find(root) == tableHeaders.end()) {
		node = node.next_sibling();
		root = node.name();
	}

#ifdef DEBUG
	std::cout << "Root is: " << root << std::endl;
	std::cout << "Now getting header keys for: " << root << std::endl;
#endif

	// Getting header keys
	const std::vector<std::string> NodeKeys = tableHeaders[root];
	Rcpp::CharacterMatrix tempRes = ret[root];

	// Prefix
	std::vector<const char*> prefix;
	prefix.resize(prefixLens[root]);

#ifdef DEBUG
	std::cout << "Applying parent" << std::endl;
#endif

	// Apply parent attributes to row and to children prefix
	for(unsigned long i = 0; i < parentPrefix.size(); i++) {
		if(parentPrefix[i] != NULL) {
			tempRes(levelCtrs[root], i) = parentPrefix[i];
			prefix[i] = parentPrefix[i];
		}
	}

	// Getting attributes
	for(pugi::xml_attribute a = node.first_attribute()
		; a
		; a = a.next_attribute()
	) {
		// Determine position
		std::string NodeKey(a.name());
		std::vector<std::string> NodeKeys = tableHeaders[root];
		unsigned col = find(NodeKeys.begin(), NodeKeys.end(), NodeKey) - NodeKeys.begin();
		if( col < NodeKeys.size() ) {
#ifdef DEBUG
			std::cout << "A:" << col << "<-" << a.value() << std::endl;

			std::cout << levelCtrs[root] << std::endl;
			std::cout << tempRes.ncol() << std::endl;
			std::cout << tempRes.nrow() << std::endl;
			std::cout << prefix.size() << std::endl;
#endif
			tempRes(levelCtrs[root], col) = a.value();
			prefix[col] = a.value();
#ifdef DEBUG
			std::cout << "Done\n" << std::endl;
#endif

		}
	}
	
	// Getting elements
	for(pugi::xml_node n = node.first_child()
		; n
		; n = n.next_sibling()
	) {
		// For echousounder's sa records
		std::string NodeKey;
		if(n.name()[0] == '\0')
			NodeKey.append(root);
		else
			NodeKey.append(n.name());

#ifdef DEBUG
		std::cout << "Transform:" << n.name() << "<->" << NodeKey << std::endl;
#endif

		// Determine position
		std::vector<std::string> NodeKeys = tableHeaders[root];
		unsigned col = find(NodeKeys.begin(), NodeKeys.end(), NodeKey) - NodeKeys.begin();

		if( col < NodeKeys.size() ) {
#ifdef DEBUG
			std::cout << "V:" << col << "<-" << n.text().as_string() << std::endl;

			std::cout << levelCtrs[root] << std::endl;
			std::cout << tempRes.ncol() << std::endl;
			std::cout << tempRes.nrow() << std::endl;
#endif
			tempRes(levelCtrs[root], col) = n.text().as_string();

			// Accomodate IDREF in the ICES XMLs
			if(strcmp(n.first_attribute().name(), "IDREF") == 0) {
				tempRes(levelCtrs[root], col) = n.first_attribute().value();
			}
#ifdef DEBUG
			std::cout << "Done\n" << std::endl;
#endif

		} else {

			// Check for missing prefix before going to children structure (for ICES XSDs)
			if(prefixLens[root] > 0 && prefix[(prefixLens[root]-1)] == NULL) {
#ifdef DEBUG
				std::cout << "We have missing prefix!" << std::endl;
				std::cout << "Peek: " << tempRes(levelCtrs[root], (prefixLens[root]-1)) << std::endl;
#endif
				for(int jj = 0; jj < prefixLens[root]; jj++) {
					if(!Rcpp::CharacterVector::is_na(tempRes(levelCtrs[root], jj)))
						prefix[jj] = tempRes(levelCtrs[root], jj);
				}
			}

			processNode(n, prefix, tableHeaders, prefixLens, levelCtrs, ret);
		}
	}

	// Increment counter
	levelCtrs[root] = levelCtrs[root] + 1;
}

// [[Rcpp::export]]
Rcpp::List readXmlCpp(Rcpp::CharacterVector inputFile, Rcpp::List xsdObjects, Rcpp::Nullable<Rcpp::CharacterVector> xsdOverride = R_NilValue, Rcpp::Nullable<Rcpp::CharacterVector> xmlEncoding = R_NilValue, bool verbose = false)
{

	pugi::xml_document doc;
	pugi::xml_node root_node;

	// Read xml using ifstream and buffer vector
	//std::ifstream iFile (inputFile[0]);
	//std::vector<char> buffer((std::istreambuf_iterator<char>(iFile)), std::istreambuf_iterator<char>());
	//buffer.push_back('\0');

	//pugi::xml_parse_result result = doc.load_buffer_inplace_own(&buffer[0], buffer.size());

	// Read file (in Windows use UTF-8 to UTF-16 conversion)
#ifndef _WIN32
        std::string filePath(inputFile[0]);
#else
        std::string filePath1(inputFile[0]);
        std::wstring filePath;
        filePath.resize(filePath1.size());
        int newSize = MultiByteToWideChar(CP_UTF8, 0, filePath1.c_str(), filePath1.length(), const_cast<wchar_t *>(filePath.c_str()), filePath1.length());
        filePath.resize(newSize);
#endif

	if (!doc.load_file(filePath.c_str())) {
		Rcpp::Rcout << "Unable to read " << inputFile[0] << std::endl;
		return -1;
	}

	// Get namespace
	char* xmlns = NULL;
	char* ns = NULL;

// We shifted all the namespace detection procedure in R
#ifdef C_DETECT_NAMESPACE
	std::string xmlStr("xmlns");
	for(pugi::xml_attribute a = doc.first_child().first_attribute()
		; a
		; a = a.next_attribute()
	) {
		std::string NodeKey(a.name());
		std::size_t found = NodeKey.find(xmlStr);
		if (found!=std::string::npos) {
			xmlns = strdup(a.value());
			// Try to get the namespace
			char *dup = strdup(NodeKey.c_str());
			strtok(dup, ":");
			char *tmp = strtok(NULL, ":");
			if(tmp != NULL)
				ns = strdup(tmp);
			free(dup);
			break;
		}
	}
#endif

	// If there is a user supplied xsd namespace
	if (xsdOverride.isNotNull()) {
		xmlns = Rcpp::as<Rcpp::CharacterVector>(xsdOverride)[0];
	}

	if (xmlns == NULL) {
		Rcpp::stop("Can not find the XML namespace, exiting...\n");
	}

	if (verbose == true) {
		Rcpp::Rcout << "Root: " <<  doc.first_child().name() << "\n";
		Rcpp::Rcout << "XML namespace: " << xmlns << "\n";
	}

	if(ns != NULL && strlen(ns) > 0 && strcmp(ns, "xsd") !=0 && strcmp(ns, "xsi") !=0) {
		Rcpp::Rcout << "XML namespace prefix: " << ns << "\n";
		Rcpp::stop("Unfortunately, namespace support is still broken!!!\n");
	} else {
		ns = NULL;
	}

	char xsd[50];
	
	// If there is a user supplied xsd namespace
	if (xsdOverride.isNotNull()) {
		sprintf (xsd, "%s.xsd", xmlns);
	} else {
		// Process namespace to get the correct XSD data
		char *token = std::strtok(xmlns, "/");

		char *one = NULL;
		char *two = token;

		while (token) {
			one = two;
			two = token;
			token = strtok(NULL, "/");
		}

		sprintf (xsd, "%s%s.xsd", one, two);
	}
	// Print out XML information
	if (verbose == true)
		Rcpp::Rcout << "Using XSD: " << xsd << std::endl;

#ifdef DEBUG
	Rcpp::Rcout << "Getting XSD objects" << std::endl;
#endif

	// Get XSD object
	Rcpp::CharacterVector root = Rcpp::as<Rcpp::List>(xsdObjects[xsd])["root"];
	Rcpp::List treeStruct = Rcpp::as<Rcpp::List>(xsdObjects[xsd])["treeStruct"];
	Rcpp::List tableHeaders = Rcpp::as<Rcpp::List>(xsdObjects[xsd])["tableHeaders"];
	Rcpp::NumericVector prefixLens = Rcpp::as<Rcpp::List>(xsdObjects[xsd])["prefixLens"];
	Rcpp::CharacterVector levelDims = Rcpp::as<Rcpp::List>(xsdObjects[xsd])["levelDims"];
	Rcpp::CharacterVector tables = Rcpp::as<Rcpp::List>(xsdObjects[xsd])["tableOrder"];

#ifdef DEBUG
	Rcpp::Rcout << "Convert headers to C++" << std::endl;
#endif

	// convert R headers to std c++
	std::vector<std::string>  tableNamesCpp;
	std::map<std::string, std::vector<std::string> > tableHeadersCpp;
	std::map<std::string, int > prefixLensCpp;

	std::string appendNS(":");
	if(ns != NULL)
		appendNS.insert(0, ns);

	for(Rcpp::CharacterVector::iterator it = tables.begin(); it != tables.end(); ++it) {
		std::string its(*it);
		std::string itsrc(*it);

		// Use Namespace for table names
		if(ns != NULL)
			its.insert(0, appendNS);

		tableNamesCpp.push_back(its);
		tableHeadersCpp[its] = Rcpp::as<std::vector<std::string> >(tableHeaders[itsrc]);

		// Appending namespace (if any) into table header names
		if(ns != NULL && tableHeadersCpp[its].size() != 0) {
			for( unsigned subit = 0; subit < tableHeadersCpp[its].size(); ++subit) {
				tableHeadersCpp[its][subit].insert(0, appendNS);
			}
		}

		prefixLensCpp[its] = prefixLens[itsrc];
	}

#ifdef DEBUG
	Rcpp::Rcout << "Finding root node" << std::endl;
#endif

	// Find our root node
	char * rootStr = root[0];
	root_node = doc.child(rootStr);

#ifdef DEBUG
	Rcpp::Rcout << "Create prefix storage" << std::endl;
#endif

	// Create prefix storage
	std::vector<const char*> prefix;

	// Allowing one level down
	std::string toErase = "";
	if(!root_node) {
		toErase = toErase + "/" + rootStr;
		Rcpp::CharacterVector downLevel = treeStruct[rootStr];
		rootStr = downLevel[0];
		root_node = doc.child(rootStr);
	}


#ifdef DEBUG
	Rcpp::Rcout << "Creating counters" << std::endl;
#endif

	// Prepare counters
	std::map<std::string, int> levelCtrs;

#ifdef DEBUG
	Rcpp::Rcout << "Creating result list and allocations" << std::endl;
#endif
	
	// Prepare the result list
	Rcpp::List ret = Rcpp::List::create();
	
	// Pre-allocations
	for(unsigned i = 0; i < tableNamesCpp.size(); i++) {
		std::string tStr(tableNamesCpp[i]);

		// Counter		
		levelCtrs[tStr.c_str()] = 0;

		// Run XPATH
		std::string te =  Rcpp::as< std::string >(levelDims[tStr.c_str()]);

		// Cut first root if needed
		if(toErase.length() > 1) {
			size_t pos = te.find(toErase);
			te.erase(pos, toErase.length());
		}

		int count = 0;

#ifdef DEBUG
		Rcpp::Rcout << "Run XPATH!" << std::endl;
#endif

		if(te != "count()") {
			pugi::xpath_query query_countnode(te.c_str());
			count = query_countnode.evaluate_number(doc);
		}

		// Matrix
		Rcpp::CharacterVector tH = tableHeaders[tStr.c_str()];

#ifdef DEBUG
		Rcpp::Rcout << te << ", "<<  count << ", " << tH.size() << std::endl;
#endif

		Rcpp::CharacterMatrix xy((int)count, tH.size());
		std::fill(xy.begin(), xy.end(), Rcpp::CharacterVector::get_na()) ;
		ret.push_back( xy, tStr );

#ifdef DEBUG
		int sz = tH.size();
		Rcpp::Rcout << "Created matrix: " << tStr.c_str() << "(" << count << ", " << sz << ", "<<  xy.size() << ")" << std::endl;
#endif
 	}

	// Naming the result list
	ret.names() = tables;

	// Process Nodes
	processNode(root_node, prefix, tableHeadersCpp, prefixLensCpp, levelCtrs, ret);

#ifdef DEBUG
	Rcpp::Rcout << "Final tally" << std::endl;
	for(int i = 0; i < tables.size(); i++) {
		std::string tStr(tables[i]);
		Rcpp::Rcout << levelCtrs[tStr] << std::endl;
	}
#endif

	// Return results and xsd name
	return Rcpp::List::create(
		Rcpp::_["xsd"]  = xsd,
		Rcpp::_["result"]  = ret
	);
}

